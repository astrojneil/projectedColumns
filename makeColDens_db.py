import yt
import numpy as np
import trident as tri
import matplotlib.pyplot as plt
from matplotlib import cm
import subprocess
from yt.data_objects.particle_filters import add_particle_filter
import h5py
from yt.units import centimeter, gram, second, Kelvin, erg
import pandas as pd
import sqlite3
import math
import io

yt.enable_parallelism()

#### Parameters ####
back = 'test_logU_neg2'
ionTable = '/Volumes/GiantDrive1/NewTridentTables/logU_neg2/logU_neg2_final.h5'
dbname = 'testlogU.db'
makeTable = False
singleFile = True
#####


#constans
kpc = 3.086e+21*centimeter
c_speed = 3.0e10  #cm/s
mp = 1.6726e-24*gram #grams
kb = 1.3806e-16*erg/Kelvin   #egs/K

#add metallicity to dataset, constant Z = 1 Zsun
def _metallicity(field, data):
    v = data['ones']  #sets metallicity to 1 Zsun
    return data.apply_units(v, "Zsun")

def adapt_array(arr):
    #http://stackoverflow.com/a/31312102/190597 (SoulNibbler)
    out = io.BytesIO()
    np.save(out, arr)
    out.seek(0)
    return sqlite3.Binary(out.read())

def convert_array(text):
    out = io.BytesIO(text)
    out.seek(0)
    return np.load(out)

#function to determine the velocity of the cloud for a particular run
#will return an array of velocities (km/s) for the appropriate velocity bins to use
#in the observational data.
def findCloudVel(directory, runName, f_list):
    velList = []
    velFrames = []
    for i in f_list:
        data = yt.load(directory+runName+'/KH_hdf5_chk_'+i)
        allDataRegion = data.all_data()
        cloudRegion = allDataRegion.cut_region(['obj["density"] >= 3.33e-25'])  #at or above 1/3 original density (1e-24)
        avg_vy_cloud = cloudRegion.quantities.weighted_average_quantity('vely', 'ones') #vely is the velocity in the radial direction (towards/away obs)

        #need to add the frame velocity!
        #get the frame velocity
        f = h5py.File(directory+runName+'/KH_hdf5_chk_'+i, 'r')
        velframe = f['real scalars'][7][1] #cm/s
        velFrames.append(velframe/1.0e5) #append frame vel in km/s
        f.close()

        vy_cloud = (avg_vy_cloud.value+velframe)/1.0e5  #convert to km/s
        velList.append(vy_cloud)

    print('Frame vels:')
    print(velFrames)
    print('Cloud vels:')
    print(velList)
    return velList, velFrames


#given one data file sphere size, and an ion, return average absorption  (given velocity in km/s)
def calcRankCol(directory, runName, runNumber, ions, velocityBin, frameVel, unRank):
    #load and add ions to the dataset
    frameVel = frameVel*1.0e5*(centimeter/second)  #convert to cm/s

    data = yt.load(directory+runName+'/KH_hdf5_chk_'+runNumber)
    data.add_field(('gas', 'metallicity'), function=_metallicity, display_name="Metallicity", units='Zsun')

    #connect to database
    conn = sqlite3.connect(dbname, detect_types=sqlite3.PARSE_DECLTYPES)
    cursor = conn.cursor()

    #add ion fields to the dataset
    for ion in ions:
        tri.add_ion_fields(data, ions=[ion['ion']], ionization_table=ionTable)

        #select entire region
        reg = data.all_data()

        #find center
        c = reg.quantities.center_of_mass()
        #c = [0.0, 0.0, 0.0]
        #set normal vector
        N = [0, 1, 0]
        #widths [right-left, top-bot, front-back] in code units
        W = [2.468e21, 2.468e21, 2.468e21]
        #set pixels along one edge
        pix = 400

        #loop through angles;
        #for y proj: r = 0
        #for x proj: r = 1
        for r in [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:

            #set normal vector
            a = math.acos(r)
            #a = d*math.pi/180.
            N = [math.cos(a), math.sin(a), 0]

            frb = yt.off_axis_projection(data, c, N, W, pix, ion['fieldname'], no_ghost=False)

            #fix issue with rotation of 0.9
            if r == 0.9:
                frb = np.rot90(frb, k =3)


            #code to plot frb if you'd like to look at it for troubleshooting purposes
            #fig, ax = plt.subplots()
            #cax = ax.matshow(np.log10(frb), interpolation='nearest', cmap=cm.get_cmap('viridis'))
            #cbar = fig.colorbar(cax)
            #fig.savefig('proj'+ion['ion']+str(r)+'.png')


            #flatten the frb to a 1d array
            flattened_num = frb.flatten()
            projectedNum = flattened_num

            projectedNum = np.array(projectedNum)

        #write the ranked column densities to a file
        #make this a column in a dataframe that you can then save as one file per run
            if singleFile:
                unRank[ion['ion']+'_'+str(r)] = projectedNum

        #write the rank column densities to the database
        #execute query to add data to table
            if makeTable:
                # Converts np.array to TEXT when inserting
                sqlite3.register_adapter(np.ndarray, adapt_array)

                # Converts TEXT to np.array when selecting
                sqlite3.register_converter("array", convert_array)
                #write to database
                if yt.is_root():
                    cursor.execute('INSERT INTO proj (coldens, ion, run_name, direction, background, timeNum) VALUES(?, ?, ?, ?, ?, ?) ', (projectedNum, ion['ion'], runName, str(r), back, velocityBin))
                    conn.commit()

    conn.close()


    return

#make tables for the database;
## !! Deletes existing tables if already defined!! ##
#fills the run table
def makeTables(runList):
    conn = sqlite3.connect(dbname)
    cursor = conn.cursor()
    cursor.execute("""DROP TABLE IF EXISTS proj;""")
    create_proj = """CREATE TABLE proj (
    id INTEGER PRIMARY KEY,
    ion TEXT,
    coldens array,
    background TEXT,
    direction TEXT,
    run_name TEXT NOT NULL,
    timeNum FLOAT,
    FOREIGN KEY (run_name) REFERENCES run (run_name)
    );"""
    cursor.execute(create_proj)


    cursor.execute("""DROP TABLE IF EXISTS run;""")
    create_run = """ CREATE TABLE run (
    run_name TEXT PRIMARY KEY,
    mach FLOAT,
    velocity FLOAT,
    tcc FLOAT
    )"""
    cursor.execute(create_run)

    for run in runList:
        exe_str = "INSERT INTO run (run_name, mach, velocity, tcc) VALUES(?, ?, ?, ?)"
        exe_param = (run['Name'], run['Mach'], run['velocity'], run['tcc'])
        cursor.execute(exe_str, exe_param)

    conn.commit()
    conn.close()

    return



def main():
##### Runs to have columns generated #######
    run1 = { 'Name':'T0.3_v1000_chi300_cond',
        'Dir':'../../Blob_paper2/Files/',
        'Mach':3.8,
        'tcc':1.7,
        'velocity':1000,
        'f_list':['0013', '0038', '0080', '0132']}
    run2 = { 'Name':'T3_v3000_chi3000_cond',
        'Dir':'../../Blob_paper2/Files/',
        'Mach':3.6,
        'tcc':1.8,
        'velocity':3000,
        'f_list':['0001', '0004', '0007', '0010']}
    run3 = { 'Name':'T1_v1700_chi1000_cond',
        'Dir':'../../Blob_paper2/Files/',
        'Mach':3.5,
        'tcc':1.8,
        'velocity':1700,
        'f_list':['0002', '0010', '0017', '0028']}
    run4 = { 'Name':'T0.3_v1000_chi300',
        'Dir':'../../Blob_paper1/Files/',
        'Mach':3.8,
        'tcc':1.7,
        'velocity':1000,
        'f_list':['0025', '0033', '0042', '0058']}
    run5 = { 'Name':'T3_v3000_chi3000',
        'Dir':'../../Blob_paper1/Files/',
        'Mach':3.6,
        'tcc':1.8,
        'velocity':3000,
        'f_list':['0021', '0030', '0040', '0062']}
    run6 = { 'Name':'T1_v1700_chi1000',
        'Dir':'../../Blob_paper1/Files/',
        'Mach':3.5,
        'tcc':1.8,
        'velocity':1700,
        'f_list':['0021', '0029', '0038', '0052']}
    run7 = { 'Name':'HC_v1000_chi300_cond',
        'Dir':'../../Blob_paper3/Files/',
        'Mach':3.8,
        'tcc':1.8,
        'velocity':1000,
        'f_list':['0054', '0060', '0080', '0107']}
    run8 = { 'Name':'HC_v1700_chi1000_cond',
        'Dir':'../../Blob_paper3/Files/',
        'Mach':3.5,
        'tcc':1.8,
        'velocity':1700,
        'f_list':['0024', '0050', '0082', '0083']}
    run9 = { 'Name':'HC_v3000_chi3000_cond',
        'Dir':'../../Blob_paper3/Files/',
        'Mach':3.6,
        'tcc':1.8,
        'velocity':3000,
        'f_list':['0007', '0015', '0026', '0049']}
    run10 = { 'Name':'LowCond_v1700_chi300_cond',
        'Dir':'../../Blob_paper3/Files/',
        'Mach':6.5,
        'tcc':1.0,
        'velocity':1700,
        'f_list':['0016', '0075', '0115', '0184']}

    run11 = { 'Name':'T0.3_v1700_chi300_cond',
        'Dir':'../../Blob_paper2/Files/',
        'Mach':6.5,
        'tcc':1.0,
        'velocity':1700,
        'f_list':['0003', '0020', '0046', '0078']}
    run12 = { 'Name':'T0.3_v3000_chi300_cond',
        'Dir':'../../Blob_paper2/Files/',
        'Mach':11.4,
        'tcc':0.56,
        'velocity':3000,
        'f_list':['0001', '0004', '0014', '0035']}
    run13 = { 'Name':'T3_v860_chi3000_cond',
        'Dir':'../../Blob_paper2/Files/',
        'Mach':1.0,
        'tcc':6.2,
        'velocity':860,
        'f_list':['0001', '0003', '0006', '0010']}
    run14 = { 'Name':'T10_v1500_chi10000_cond',
        'Dir':'../../Blob_paper2/Files/',
        'Mach':1.0,
        'tcc':6.5,
        'velocity':1500,
        'f_list':['0001', '0002', '0004', '0008']}
    run15 = { 'Name':'T1_v480_chi1000_cond',
        'Dir':'../../Blob_paper2/Files/',
        'Mach':1.0,
        'tcc':6.4,
        'velocity':480,
        'f_list':['0002', '0009', '0017', '0031']}

    run16 = { 'Name':'T0.3_v1700_chi300',
        'Dir':'../../Blob_paper1/Files/',
        'Mach':6.5,
        'tcc':1.0,
        'velocity':1700,
        'f_list':['0022', '0032', '0053', '0085']}
    run17 = { 'Name':'T0.3_v3000_chi300',
        'Dir':'../../Blob_paper1/Files/',
        'Mach':11.4,
        'tcc':0.56,
        'velocity':3000,
        'f_list':['0028', '0044', '0065', '0110']}
    run18 = { 'Name':'T3_v430_chi3000',
        'Dir':'../../Blob_paper1/Files/',
        'Mach':0.5,
        'velocity':430,
        'tcc':12.5,
        'f_list':['0010', '0011', '0016', '0024']}
    run19 = { 'Name':'T3_v860_chi3000',
        'Dir':'../../Blob_paper1/Files/',
        'Mach':1.0,
        'velocity':860,
        'tcc':6.2,
        'f_list':['0010', '0018', '0030', '0038']}
    run20 = { 'Name':'T1_v3000_chi1000',
        'Dir':'../../Blob_paper1/Files/',
        'Mach':6.2,
        'velocity':3000,
        'tcc':1.0,
        'f_list':['0022', '0032', '0048', '0095']}
    run21 = { 'Name':'T10_v1500_chi10000',
        'Dir':'../../Blob_paper1/Files/',
        'Mach':1.0,
        'velocity':1500,
        'tcc':6.5,
        'f_list':['0014', '0021', '0029', '0044']}
    run22 = { 'Name':'T1_v480_chi1000',
        'Dir':'../../Blob_paper1/Files/',
        'Mach':1.0,
        'velocity':480,
        'tcc':6.4,
        'f_list':['0009', '0013', '0024', '0035']}
    run23 = { 'Name':'T1_v1700_chi1000_lref6',
        'Dir':'../../Blob_paper1/Files/',
        'Mach':3.5,
        'velocity':1700,
        'tcc':1.8,
        'f_list':['0021', '0029', '0038']}
    run24 = { 'Name':'CoolFloor2e4',
        'Dir':'/Volumes/GiantDrive2/',
        'Mach':6.5,
        'velocity':1700,
        'tcc':1.0,
        'f_list':['0010', '0020', '0098', '0103']}


#add the runs to the list that will have columns ranked
    runList = []
    #add only the 3 different conduction levels
    runList.append(run10)
    runList.append(run11)
    runList.append(run16)


    #other potential runs to use:
    '''
    run6 = { 'Name':'T1_v1700_chi1000_lref4',
                'Dir':'../../Blob_paper1/Files/',
                'Mach':3.5,
                'f_list':['0020', '0030', '0040', '0050']}
    run7 = { 'Name':'T1_v1700_chi1000_lref6',
                'Dir':'../../Blob_paper1/Files/',
                'Mach':3.5,
                'f_list':['0021', '0029', '0038']}
    '''


### dictionaries of ion info
    ion1 = {'ion':'O VI',
        'fieldname':'O_p5_number_density',
        'ionfolder': '/OVI/',
        'rest_wave': 1031.91,
        'sigma': 1.1776e-18,
        'massNum': 16.0}
    ion2 = {'ion':'C IV',
        'fieldname':'C_p3_number_density',
        'ionfolder': '/CIV/',
        'rest_wave': 1548.18,
        'sigma': 2.5347e-18,
        'massNum': 12.0}
    ion3 = {'ion':'N V',
        'fieldname':'N_p4_number_density',
        'ionfolder': '/NV/',
        'rest_wave': 1242.8,
        'sigma': 8.3181e-19,
        'massNum': 14.0}
    ion4 = {'ion':'C II',
        'fieldname':'C_p1_number_density',
        'ionfolder': '/CII/',
        'rest_wave': 1335.66,
        'sigma': 1.4555e-19,
        'massNum': 12.0}
    ion5 = {'ion': 'Ne VIII',
            'fieldname': 'Ne_p7_number_density',
            'ionfolder':'/NeVIII/',
            'rest_wave': 770.406,
            'sigma': 6.74298e-19,
            'massNum': 20.0}
    ion6 = {'ion': 'C III',
            'fieldname': 'C_p2_number_density',
            'ionfolder':'/CIII/',
            'rest_wave': 977.02,
            'sigma': 6.359e-18,
            'massNum': 12.0}
    ion7 = {'ion': 'Mg II',
            'fieldname': 'Mg_p1_number_density',
            'ionfolder':'/MgII/',
            'rest_wave': 1239.92,
            'sigma': 6.60717e-21,
            'massNum': 24.0}
    ion8 = {'ion': 'Si III',
            'fieldname': 'Si_p2_number_density',
            'ionfolder':'/SiIII/',
            'rest_wave': 1206.5,
            'sigma': 2.2258e-20,
            'massNum': 28.0}
    ion9 = {'ion': 'Si IV',
            'fieldname': 'Si_p3_number_density',
            'ionfolder':'/SiIV/',
            'rest_wave': 1402.77,
            'sigma': 3.0694e-18,
            'massNum': 28.0}
    ion10 = {'ion': 'H I 1215',
            'fieldname': 'H_p0_number_density',
            'ionfolder':'/HI/',
            'rest_wave': 1215.67,
            'sigma': 4.3394e-18,
            'massNum': 2.0}
    ion11 = {'ion': 'H II',
            'fieldname': 'H_p1_number_density',
            'ionfolder':'/HII/',
            'rest_wave': 1215.67,
            #'sigma': 4.3394e-18,
            'massNum': 2.0}

#create list of IOns
    ionList = []
    ionList.append(ion1)
    ionList.append(ion10)
    ionList.append(ion2)
    ionList.append(ion3)
    ionList.append(ion4)
    ionList.append(ion5)
    ionList.append(ion6)
    ionList.append(ion7)
    ionList.append(ion8)
    ionList.append(ion9)
    ionList.append(ion11)

    #make database if you want to make tables
    if makeTables:
        makeTables(runList)


### run though functions!
    for run in runList:
        #dataframe for a single file
        unRank = pd.DataFrame()

        #find appropriate velocity bins
        velBins, velFrames = findCloudVel(run['Dir'], run['Name'], run['f_list'])

      #velocities are in km/s
        for v in range(len(velBins)): #velocity = velBins[v]  for t50: [velBins[2]]
            #calculate the ranked columns for each ion/velocity bin
            calcRankCol(run['Dir'], run['Name'], run['f_list'][v], ionList, v, velFrames[v], unRank)

        #save a single file for the run
        if singleFile:
            unRank.to_csv(run['Name']+'_y_HM1e7.csv')

if __name__ =="__main__":
    main()
