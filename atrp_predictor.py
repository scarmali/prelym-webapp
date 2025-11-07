#!/usr/bin/env python
# coding: utf-8

# In[1]:


from biopandas.pdb import PandasPdb
import pandas as pd
from Bio.PDB import *
import freesasa
from MDAnalysis import *
# from propkatraj import get_propka  # Not available in current version
import itertools
import numpy as np
# Robust DSSP function with fallbacks
def dsspf(filename):
    """DSSP function with multiple fallback strategies"""
    try:
        # First try ssbio's DSSP
        from ssbio.protein.structure.properties.dssp import get_dssp_df_on_file
        result = get_dssp_df_on_file(filename)

        # Check if result is empty or missing required columns
        if result.empty or 'resnum' not in result.columns:
            raise Exception("DSSP returned empty or invalid result")

        return result

    except Exception as e:
        print(f"Warning: ssbio DSSP failed ({e}), trying BioPython DSSP...")

        try:
            # Try BioPython's DSSP
            from Bio.PDB import PDBParser
            from Bio.PDB.DSSP import DSSP
            parser = PDBParser()
            structure = parser.get_structure("protein", filename)
            model = structure[0]

            dssp = DSSP(model, filename)
            data = []
            for key in dssp.keys():
                chain_id = key[0]
                res_id = key[1][1]
                sec_struct = dssp[key][2]
                # Keep original DSSP codes in 'aa' column and add interpreted structure
                if sec_struct in ['H', 'G', 'I']:
                    ss = 'Helix'
                elif sec_struct in ['B', 'E']:
                    ss = 'Strand'
                else:
                    ss = 'Coil'
                data.append({'resnum': res_id, 'sec_struc': ss, 'chain': chain_id, 'aa': sec_struct})
            return pd.DataFrame(data)

        except Exception as e2:
            print(f"Warning: BioPython DSSP also failed ({e2}), using dihedral-based prediction...")

            # Final fallback: dihedral angle-based secondary structure prediction
            try:
                from Bio.PDB import PDBParser, calc_dihedral
                import numpy as np
                parser = PDBParser()
                structure = parser.get_structure("protein", filename)

                data = []
                for model in structure:
                    for chain in model:
                        chain_id = chain.get_id()
                        residues = [res for res in chain if res.get_id()[0] == ' ']

                        for i, residue in enumerate(residues):
                            res_num = residue.get_id()[1]
                            ss = 'Coil'  # Default
                            aa_code = 'C'

                            try:
                                # Calculate phi and psi angles for secondary structure prediction
                                phi = None
                                psi = None

                                if i > 0 and i < len(residues) - 1:
                                    prev_res = residues[i-1]
                                    next_res = residues[i+1]

                                    # Calculate phi angle (C(i-1) - N(i) - CA(i) - C(i))
                                    if ('C' in prev_res and 'N' in residue and
                                        'CA' in residue and 'C' in residue):
                                        phi = calc_dihedral(
                                            prev_res['C'].get_vector(),
                                            residue['N'].get_vector(),
                                            residue['CA'].get_vector(),
                                            residue['C'].get_vector()
                                        )
                                        phi = np.degrees(phi)

                                    # Calculate psi angle (N(i) - CA(i) - C(i) - N(i+1))
                                    if ('N' in residue and 'CA' in residue and
                                        'C' in residue and 'N' in next_res):
                                        psi = calc_dihedral(
                                            residue['N'].get_vector(),
                                            residue['CA'].get_vector(),
                                            residue['C'].get_vector(),
                                            next_res['N'].get_vector()
                                        )
                                        psi = np.degrees(psi)

                                # Improved DSSP-like Ramachandran classification
                                if phi is not None and psi is not None:
                                    # Alpha helix region (more precise boundaries)
                                    if (-80 <= phi <= -40) and (-60 <= psi <= -20):
                                        ss = 'Helix'
                                        aa_code = 'H'
                                    # 310-helix region (tighter helix)
                                    elif (-90 <= phi <= -40) and (-45 <= psi <= 0):
                                        ss = 'Helix'
                                        aa_code = 'G'
                                    # Beta sheet regions (antiparallel and parallel)
                                    elif (-180 <= phi <= -90) and (90 <= psi <= 180):
                                        ss = 'Strand'
                                        aa_code = 'E'
                                    elif (-180 <= phi <= -90) and (-180 <= psi <= -100):
                                        ss = 'Strand'
                                        aa_code = 'E'
                                    # Extended beta region
                                    elif (-150 <= phi <= -90) and (100 <= psi <= 180):
                                        ss = 'Strand'
                                        aa_code = 'E'
                                    # Polyproline II and left-handed helix regions
                                    elif (-90 <= phi <= -40) and (120 <= psi <= 180):
                                        ss = 'Coil'
                                        aa_code = 'C'
                                    # Beta turn regions
                                    elif (-90 <= phi <= 0) and (-40 <= psi <= 40):
                                        ss = 'Turn'
                                        aa_code = 'T'
                                    # Everything else is coil
                                    else:
                                        ss = 'Coil'
                                        aa_code = 'C'
                            except:
                                # If angle calculation fails, default to coil
                                ss = 'Coil'
                                aa_code = 'C'

                            data.append({'resnum': res_num, 'sec_struc': ss, 'chain': chain_id, 'aa': aa_code})

                return pd.DataFrame(data) if data else pd.DataFrame({'resnum': [1], 'sec_struc': ['Coil'], 'chain': ['A'], 'aa': ['C']})

            except Exception as e3:
                print(f"Warning: All DSSP methods failed ({e3}), returning minimal structure...")
                # Last resort: return minimal structure
                return pd.DataFrame({'resnum': [1], 'sec_struc': ['Coil'], 'chain': ['A'], 'aa': ['C']})
from biopandas.pdb import PandasPdb
from propkatraj import PropkaTraj
import MDAnalysis as mda
import mdtraj as md
import math
from collections import defaultdict


# In[2]:


def pdb2esa(filename,probe_radius):

    #ESA calculation:
    parser = PDBParser()
    structure = parser.get_structure('self',filename)
    param = freesasa.Parameters()
    param.setProbeRadius(probe_radius)
    result, sasa_classes = freesasa.calcBioPDB(structure,param)

    print(f'The defined probe radius is {param.probeRadius()} Å.')

    # Determine how many atoms FreeSASA actually processed
    freesasa_atom_count = 0
    try:
        # Find the maximum valid atom index
        for i in range(10000):  # reasonable upper limit
            try:
                result.atomArea(i)
                freesasa_atom_count = i + 1
            except:
                break
    except:
        freesasa_atom_count = 0

    print(f'FreeSASA processed {freesasa_atom_count} atoms')

    #Biopython:
    ppdb = PandasPdb().read_pdb(filename)
    df = ppdb.df
    df = df['ATOM']

    # Only use the atoms that FreeSASA processed
    df = df.head(freesasa_atom_count)

    atom_no = df['atom_number']
    chain = df['chain_id']
    resid_no = df['residue_number']
    resid_no = resid_no.astype(str)

    chain_list = []

    for i in range(len(chain)):
        val = resid_no.iloc[i]+'.'+chain.iloc[i]
        chain_list.append(val)

    df = pd.DataFrame({'Chain':chain_list,'Atom Number':atom_no})
    residue_names = df['Chain'].unique()

    esa = []
    a = 0

    for i in range(len(residue_names)):
        total = 0
        df_new = df[df['Chain'] == residue_names[i]]
        length = len(df_new)

        for j in range(a,a+length):
            if j < freesasa_atom_count:
                area = result.atomArea(j)
                total += area
            else:
                print(f'Warning: Skipping atom {j} as it exceeds FreeSASA processed atoms')

        a += length
        esa.append(total)
        
    ESA = pd.DataFrame({'Chain':residue_names,'ESA':esa})
    
    residue = []
    chain = []
    esa_val = []

    for i in range(len(ESA)):
        res,cha = ESA['Chain'][i].split('.') 
        val = ESA['ESA'][i]
        residue.append(res)
        chain.append(cha)
        esa_val.append(val)
    data = pd.DataFrame({'Residue':residue,'Chain':chain,'ESA':esa_val})
    data['Residue'] = data['Residue'].astype(int)
    
    return data


# In[3]:


def pdb2pka(filename):
    
    u = mda.Universe(filename)
    pkatraj = PropkaTraj(u, select='protein', skip_failure=False)
    pkatraj.run()
    df = pkatraj.results.pkas
    index = df.T[0].index.values
    values = df.T[0].values
    df = pd.DataFrame({'Residue':index,'pKa':values})
    
    return df


# In[4]:


def pdb2hbond(filename,filename1):
    
    t = md.load_pdb(filename1)
    hbonds = md.baker_hubbard(t, periodic=False)
    label = lambda hbond : '%s -- %s' % (t.topology.atom(hbond[0]), t.topology.atom(hbond[2]))
    
    hbonds_list = []

    for hbond in hbonds:
        val = label(hbond)
        hbonds_list.append(val)
        
    resid_hbonds = []

    for i in range(len(hbonds_list)):
        a = hbonds_list[i].split('-')[0]
        p,q = a[:3],a[3:]
    
        if p == 'LYS':
            resid_hbonds.append(q)
            
    df = pd.DataFrame({'Residue':resid_hbonds})
    df['Residue'] = df['Residue'].astype(int)
    
    ppdb = PandasPdb().read_pdb(filename)
    df1 = ppdb.df
    df2 = df1['ATOM']
    chains = df2['chain_id'].unique()
    min_chain = []
    max_chain = []
    
    for i in range(len(chains)):
        mini = df2[df2['chain_id'] == chains[i]]['residue_number'].min()
        maxi = df2[df2['chain_id'] == chains[i]]['residue_number'].max()
        min_chain.append(mini)
        max_chain.append(maxi)
        
    index = np.zeros(len(df))
    counter = 0

    for i in range(1,len(df)-1):
        if df['Residue'][i] < df['Residue'][i-1]:
            counter += 1
        index[i] = counter
    
    index[-1] = counter
    
    df['Counter'] = index
    df['Counter'] = df['Counter'].astype(int)
    spl = int(len(chains)/len(df['Counter'].unique()))
    
    chain_name = []
    count = 0
    
    for i in range(len(df['Counter'].unique())):
        data = df[df['Counter'] == i]
        data = data.reset_index(drop=True)
        name = []
        
        for j in range(len(data)):
            for k in range(count,count+spl):
                if (data['Residue'][j]>min_chain[k]) and (data['Residue'][j]<max_chain[k]):
                    name = chains[k]
            
            chain_name.append(name)
        count += spl
        
    df['Chain'] = chain_name
    df = df.drop(['Counter'],axis=1)
    
    return df


# In[5]:


def pdb2ss(filename):
    
    ppdb = PandasPdb().read_pdb(filename)
    df = ppdb.df
    df = df['ATOM']
    df = df[df['residue_name'] == 'LYS']
    df = df.drop_duplicates(subset=['chain_id','residue_number'])
    df = df.reset_index(drop=True)
    a = df['chain_id']
    b = df['residue_number']
    resid_data = pd.DataFrame({'Chain':a,'Residue':b})
    ss_data = dsspf(filename)
    ss_data = ss_data.reset_index(drop=True)
    
    ss_list = []

    for i in range(len(resid_data)):
    
        if ss_data[ss_data['resnum'] == resid_data['Residue'][i]].empty:
            ss__ = 'Coil'
        
        else:
            for j in range(len(ss_data)):
                
                if (resid_data['Residue'][i] == ss_data['resnum'][j]) and (resid_data['Chain'][i]==ss_data['chain'][j]):
            
                    ss_ = ss_data['aa'][j]
            
                    if ss_ == 'H':
                        ss__ = 'Helix'
                    
                    elif ss_ == 'B':
                        ss__ = 'Coil'
                    
                    elif ss_ == 'E':
                        ss__ = 'Strand'
                    
                    elif ss_ == 'G':
                        ss__ = 'Helix'
                    
                    elif ss_ == 'I':
                        ss__ = 'Helix'
                    
                    elif ss_ == 'T':
                        ss__ = 'Coil'
                    
                    elif ss_ == 'S':
                        ss__ = 'Coil'
                
                    else :
                        ss__ = 'Coil'
                
                
        ss_list.append(ss__)
        
    resid_data['Secondary Structure'] = ss_list
    
    return resid_data


# In[6]:


def pdb2res(filename):
    
    ppdb = PandasPdb().read_pdb(filename)
    df = ppdb.df
    data = df['ATOM']
    data = data.drop_duplicates(subset=['chain_id','residue_number'])
    data = data.reset_index(drop=True)
    residue = []
    chain = []
    amino_group = []

    for i in range(len(data)):
        res = data['residue_number'][i]
        cha = data['chain_id'][i]
        amin = data['residue_name'][i]
    
        residue.append(res)
        chain.append(cha)
        amino_group.append(amin)

    df = pd.DataFrame({'Chain':chain,'Residue':residue,'Amino Group':amino_group})
    n_termi = df.drop_duplicates('Chain')
    N = ['N'] * len(n_termi)
    n_termi = n_termi.drop(['Amino Group'],axis=1)
    n_termi['Amino Group'] = N
    lys_ = df[df['Amino Group'] == 'LYS']
    df = pd.concat([n_termi, lys_], ignore_index=True)
    df = df.sort_values(by=['Chain','Residue'])
    df = df.reset_index(drop=True)
    
    return df


# In[7]:


def pka_final(filename):
    
    data = pdb2res(filename)
    num = data.groupby('Residue').count().min().min()
    pka_data = pdb2pka(filename)
    min_max = data.groupby('Chain').agg({'Residue':['min','max']}).reset_index()
    
    index_list = defaultdict(list)
    cur_min = min_max.iloc[0,1]
    cur_max = min_max.iloc[0,2]
    i=0
    isinrange = True

    for index,row in pka_data.iterrows():

        if (row['Residue'] >= cur_min) and (row['Residue'] <= cur_max) and isinrange:

            index_list[min_max.iloc[i,0]].append(index)

            if row['Residue'] == cur_max:
                isinrange = False

        elif min_max.shape[0]-1 > i:

            if (row['Residue'] >= min_max.iloc[i+1,1]) and (row['Residue'] <= min_max.iloc[i+1,2]):

                i += 1   
                cur_min = min_max.iloc[i,1]
                cur_max = min_max.iloc[i,2]

                index_list[min_max.iloc[i,0]].append(index)

                isinrange = True

            else:
                isinrange = False

        else:
            continue
            
    pka_data['Chain'] = ''
    for chain,index in index_list.items():
        pka_data.iloc[index,2] = chain
        
    data = pka_data
    data2 = pdb2res(filename)
    data3 = pd.concat([data, data2])
    data3 = data3.drop(columns=['Amino Group'])
    data3 = data3.sort_values(['Chain','Residue'])
    index_vals = data3[data3['Chain'] == ''].index
    data3 = data3.drop(index_vals).reset_index(drop=True)
    data3['Amino Group'] = 'LYS'
    index_ = list(data3.drop_duplicates('Chain').index)
    data3.loc[index_,'Amino Group'] = 'N'
    data3 = data3.drop_duplicates(['Chain','Residue','Amino Group'])
    data3 = data3.reset_index(drop=True)

    data2['pKa'] = 0

    for i in range(len(data2)):
        for j in range(len(data3)):

            if (data2['Chain'][i]==data3['Chain'][j]) and (data2['Residue'][i]==data3['Residue'][j]) and (data2['Amino Group'][i]==data3['Amino Group'][j]):
                data2.loc[i,'pKa'] = data3['pKa'][j]

    data2['pKa'] = data2['pKa'].replace(np.nan,0)
    
    return data2


# In[8]:


def search(filename,list_strings):
    
    line_number = 0
    results = []
    
    with open(filename,'r') as read_obj:
        
        for line in read_obj:
            line_number += 1
            
            for string_to_search in list_strings:
                if string_to_search in line:
                    results.append(line_number)
                             
    return results


# In[9]:


def pdb2charge(filename,filename2):
    
    line_no = search(filename2,['ATOM'])
    max_line = np.max(line_no)
    min_line = np.min(line_no)
    data = pd.read_csv(filename2,skiprows=min_line-1,nrows = max_line-min_line,delimiter=r'\s+',header=None)
    data = data.rename(columns={1:'Atom No.',2:'Atom',3:'Amino Group',4:"Residue",5:'X',6:'Y',7:'Z',8:'Charge',9:'Radius'})

    lys_data = data[data['Amino Group'] == 'LYS']
    lys_data = lys_data[lys_data['Atom']=='NZ']

    surface_charge = np.zeros(len(data))

    for i in lys_data.index:
    
        value = 0
        X1 = data['X'][i]
        Y1 = data['Y'][i]
        Z1 = data['Z'][i]
        q1 = data['Charge'][i]
    
        for j in range(len(data)):
        
            X2 = data['X'][j]
            Y2 = data['Y'][j]
            Z2 = data['Z'][j]
            q2 = data['Charge'][j]
            r = np.sqrt((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2)
    
            if r > 0:
                F = -333*q2*q1/(r)
                value += F
            
        surface_charge[i] = value    
    
    data['Surface Charge'] = surface_charge
    data = data[data['Amino Group'] == 'LYS']
    data = data[data['Atom'] == 'NZ']
    data = data.reset_index(drop=True)
    local_charge = [''] * len(data)

    for i in range(len(data)):
    
        if data['Surface Charge'][i] >100:
            local_charge[i] = 'Yes'
    
        else:
            local_charge[i] = 'No'
        
    data['Area Of Lower Charge'] = local_charge
    
    res = pdb2res(filename)
    charge1 = data.drop_duplicates('Residue',keep='first')
    charge2 = data.drop_duplicates('Residue',keep='last')
    res1 = res.drop_duplicates('Residue',keep='first')
    res2 = res.drop_duplicates('Residue',keep='last')
    res1 = res1.reset_index(drop=True)
    res2 = res2.reset_index(drop=True)
    charge1 = charge1.reset_index(drop=True)
    charge2 = charge2.reset_index(drop=True)

    cha_1 = ['']*len(res1)
    cha_2 = ['']*len(res2)

    for i in range(len(res1)):
        for j in range(len(charge1)):
        
            if res1['Residue'][i] == charge1['Residue'][j]:
                cha_1[i] = charge1['Area Of Lower Charge'][j]
            
    for i in range(len(res2)):
        for j in range(len(charge2)):
        
            if res2['Residue'][i] == charge2['Residue'][j]:
                cha_2[i] = charge2['Area Of Lower Charge'][j]
            
    res1['Area Of Lower Charge'] = cha_1
    res2['Area Of Lower Charge'] = cha_2
    data = pd.concat([res1, res2], ignore_index=True)
    data = data.reset_index(drop=True)
    data['Area Of Lower Charge'] = data['Area Of Lower Charge'].replace('',np.nan)
    data
    
    return data


# In[10]:


def final_data_table(filename,filename1,filename2,probe_radius):
    
    resids = pdb2res(filename)
    esa_data = pdb2esa(filename,probe_radius)

    esa_vals = np.zeros(len(resids))

    for i in range(len(resids)):
        for j in range(len(esa_data)):
        
            if (resids['Residue'][i] == esa_data['Residue'][j]) and (resids['Chain'][i] == esa_data['Chain'][j]):
                esa_vals[i] = esa_data['ESA'][j]

    main_data = resids
    main_data['ESA'] = esa_vals
    pka_vals = pka_final(filename)
    pkavals = np.zeros(len(main_data))

    for i in range(len(main_data)):
        for j in range(len(pka_vals)):
        
            if (main_data['Residue'][i]==pka_vals['Residue'][j]) and (main_data['Chain'][i]==pka_vals['Chain'][j]) and (main_data['Amino Group'][i]==pka_vals['Amino Group'][j]):
                pkavals[i] = pka_vals['pKa'][j]
            
    main_data['pKa'] = pkavals
    main_data['pKa'] = main_data['pKa'].replace(0.0,np.nan)

    ss_data = pdb2ss(filename)
    ss_vals = ['']*len(main_data)

    for i in range(len(main_data)):
        for j in range(len(ss_data)):
        
            if (main_data['Residue'][i] == ss_data['Residue'][j]) and (main_data['Chain'][i] == ss_data['Chain'][j]):
                ss_vals[i] = ss_data['Secondary Structure'][j]
            
    main_data['Secondary Structure'] = ss_vals
    main_data['Secondary Structure'] = main_data['Secondary Structure'].replace('',np.nan)
    hbond_data = pdb2hbond(filename,filename1)
    hbond_vals = ['']*len(main_data)

    for i in range(len(main_data)):
        for j in range(len(hbond_data)):
        
            if (main_data['Residue'][i] == hbond_data['Residue'][j]) and (main_data['Chain'][i] == hbond_data['Chain'][j]):
                hbond_vals[i] = 'Yes'
            
    main_data['H-Donor'] = hbond_vals
    main_data['H-Donor'] = main_data['H-Donor'].replace('','No')

    charge_data = pdb2charge(filename,filename2)

    charge_val = ['']*len(main_data)

    for i in range(len(main_data)):
        for j in range(len(charge_data)):
        
            if (main_data['Residue'][i]==charge_data['Residue'][j]) and (main_data['Amino Group'][i]==charge_data['Amino Group'][j]):
                charge_val[i] = charge_data['Area Of Lower Charge'][j]
            
    main_data['Area Of Lower Charge'] = charge_val
    
    return main_data


# In[11]:


def decision_tree(filename,filename1,filename2,probe_radius):
    
    df = final_data_table(filename,filename1,filename2,probe_radius)

    interaction = []

    chain_names = df.Chain.unique()

    df_chain = []

    for name in chain_names:
    
        df_ = df[df['Chain'] == name]
        df_ = df_.reset_index(drop=True)
        df_chain.append(df_)

    for i in range(len(df_chain)):
    
        for j in range(len(df_chain[i])):
        
            #ESA filter

            if df_chain[i]['ESA'][j] <= 50:
                inter = 'non-reacting'

            elif df_chain[i]['ESA'][j] > 50:
            
                #First Entry:
                
                if j == 0:
                
                    if len(df_chain[i]) < 2:
                    
                        if df_chain[i]['Amino Group'][j] == 'N':
                            inter = 'fast-reacting'
                                
                        else:
                            
                            if df_chain[i]['Secondary Structure'][j] == 'Helix':
                                inter = 'slow-reacting'
                            
                            elif df_chain[i]['Secondary Structure'][j] == 'Strand':
                                if df_chain[i]['pKa'][j] <= 10.3:
                                    if df_chain[i]['Area Of Lower Charge'][j] == 'Yes':
                                        inter = 'fast-reacting'
                                    else:
                                        inter = 'slow-reacting'
                                    
                                else:
                                    inter = 'slow-reacting'
                                
                            else:  # Coil, Turn, or other flexible structures
                                if df_chain[i]['ESA'][j] >= 100:
                                    if df_chain[i]['Area Of Lower Charge'][j] == 'Yes':
                                        inter = 'fast-reacting'
                                    
                                    elif df_chain[i]['H-Donor'][j] == 'Yes':
                                        inter = 'fast-reacting'
                                    
                                    else:
                                        inter = 'slow-reacting'
                                    
                                else:
                                    inter = 'slow-reacting'
                

                    
                
                    #Steric Hindrance:
                
                    elif (df_chain[i]['Residue'][j] == df_chain[i]['Residue'][j+1] - 1) and (df_chain[i]['ESA'][j+1] > 50):
                        if df_chain[i]['pKa'][j] > df_chain[i]['pKa'][j+1]:
                            inter = 'non-reacting'
                        
                        else:
                        
                            if df_chain[i]['Amino Group'][j] == 'N':
                                inter = 'fast-reacting'
                            
                            else:
                            
                                if df_chain[i]['Secondary Structure'][j] == 'Helix':
                                    inter = 'slow-reacting'
                                
                                elif df_chain[i]['Secondary Structure'][j] == 'Strand':
                                    if df_chain[i]['pKa'][j] <= 10.3:
                                        if df_chain[i]['Area Of Lower Charge'][j] == 'Yes':
                                            inter = 'fast-reacting'
                                        else:
                                            inter = 'slow-reacting'
                                    
                                    else:
                                        inter = 'slow-reacting'
                                
                                else:
                                    if df_chain[i]['ESA'][j] >= 100:
                                        if df_chain[i]['Area Of Lower Charge'][j] == 'Yes':
                                            inter = 'fast-reacting'
                                    
                                        elif df_chain[i]['H-Donor'][j] == 'Yes':
                                            inter = 'fast-reacting'
                                    
                                        else:
                                            inter = 'slow-reacting'
                                        
                                    else:
                                        inter = 'slow-reacting'
                                            
                                            
                                            
                    else:
                    
                        if df_chain[i]['Amino Group'][j] == 'N':
                            inter = 'fast-reacting'
                            
                        else:
                        
                            if df_chain[i]['Secondary Structure'][j] == 'Helix':
                                inter = 'slow-reacting'
                                
                            elif df_chain[i]['Secondary Structure'][j] == 'Strand':
                                if df_chain[i]['pKa'][j] <= 10.3:
                                    if df_chain[i]['Area Of Lower Charge'][j] == 'Yes':
                                        inter = 'fast-reacting'
                                    else:
                                        inter = 'slow-reacting'
                                    
                                else:
                                    inter = 'slow-reacting'
                                
                            else:  # Coil, Turn, or other flexible structures
                                if df_chain[i]['ESA'][j] >= 100:
                                    if df_chain[i]['Area Of Lower Charge'][j] == 'Yes':
                                        inter = 'fast-reacting'
                                    
                                    elif df_chain[i]['H-Donor'][j] == 'Yes':
                                        inter = 'fast-reacting'
                                    
                                    else:
                                        inter = 'slow-reacting'
                                    
                                else:
                                    inter = 'slow-reacting'
                                        
                                        
                #Last Entry:
            
                elif j == len(df_chain[i]) - 1 :
                
                    #Steric Hindrance:
                
                    if (df_chain[i]['Residue'][j] == df_chain[i]['Residue'][j-1] + 1) and (df_chain[i]['ESA'][j-1] > 50):
                        if df_chain[i]['pKa'][j] > df_chain[i]['pKa'][j-1]:
                            inter = 'non-reacting'
                        
                        else:
                        
                            if df_chain[i]['Amino Group'][j] == 'N':
                                inter = 'fast-reacting'
                            
                            else:
                            
                                if df_chain[i]['Secondary Structure'][j] == 'Helix':
                                    inter = 'slow-reacting'
                                
                                elif df_chain[i]['Secondary Structure'][j] == 'Strand':
                                    if df_chain[i]['pKa'][j] <= 10.3:
                                        if df_chain[i]['Area Of Lower Charge'][j] == 'Yes':
                                            inter = 'fast-reacting'
                                        else:
                                            inter = 'slow-reacting'
                                    
                                    else:
                                        inter = 'slow-reacting'
                                
                                else:
                                    if df_chain[i]['ESA'][j] >= 100:
                                        if df_chain[i]['Area Of Lower Charge'][j] == 'Yes':
                                            inter = 'fast-reacting'
                                    
                                        elif df_chain[i]['H-Donor'][j] == 'Yes':
                                            inter = 'fast-reacting'
                                    
                                        else:
                                            inter = 'slow-reacting'
                                        
                                    else:
                                        inter = 'slow-reacting'
                                            
                                                                      
                    else:
                    
                        if df_chain[i]['Amino Group'][j] == 'N':
                            inter = 'fast-reacting'
                            
                        else:
                        
                            if df_chain[i]['Secondary Structure'][j] == 'Helix':
                                inter = 'slow-reacting'
                                
                            elif df_chain[i]['Secondary Structure'][j] == 'Strand':
                                if df_chain[i]['pKa'][j] <= 10.3:
                                    if df_chain[i]['Area Of Lower Charge'][j] == 'Yes':
                                        inter = 'fast-reacting'
                                    else:
                                        inter = 'slow-reacting'
                                    
                                else:
                                    inter = 'slow-reacting'
                                    
                            else:  # Coil, Turn, or other flexible structures
                                if df_chain[i]['ESA'][j] >= 100:
                                    if df_chain[i]['Area Of Lower Charge'][j] == 'Yes':
                                        inter = 'fast-reacting'
                                    
                                    elif df_chain[i]['H-Donor'][j] == 'Yes':
                                        inter = 'fast-reacting'
                                    
                                    else:
                                        inter = 'slow-reacting'
                                    
                                else:
                                    inter = 'slow-reacting'
                                        


                #Middle Entries:
            
                elif (j>0) and (j<len(df_chain[i])-1):
                
                    #Previous Entry:
                
                    if (df_chain[i]['Residue'][j] == df_chain[i]['Residue'][j-1] + 1) and (df_chain[i]['ESA'][j-1]> 50):
                    
                        if df_chain[i]['pKa'][j] > df_chain[i]['pKa'][j-1]:
                            inter = 'non-reacting'
                        
                        else:
                        
                            if df_chain[i]['Amino Group'][j] == 'N':
                                inter = 'fast-reacting'
                            
                            else:
                            
                                if df_chain[i]['Secondary Structure'][j] == 'Helix':
                                    inter = 'slow-reacting'
                                
                                elif df_chain[i]['Secondary Structure'][j] == 'Strand':
                                    if df_chain[i]['pKa'][j] <= 10.3:
                                        if df_chain[i]['Area Of Lower Charge'][j] == 'Yes':
                                            inter = 'fast-reacting'
                                        else:
                                            inter = 'slow-reacting'
                                    
                                    else:
                                        inter = 'slow-reacting'
                                
                                else:
                                    if df_chain[i]['ESA'][j] >= 100:
                                        if df_chain[i]['Area Of Lower Charge'][j] == 'Yes':
                                            inter = 'fast-reacting'
                                    
                                        elif df_chain[i]['H-Donor'][j] == 'Yes':
                                            inter = 'fast-reacting'
                                    
                                        else:
                                            inter = 'slow-reacting'
                                        
                                    else:
                                        inter = 'slow-reacting'
                                            
                
                    #Later Entry:
                
                    elif (df_chain[i]['Residue'][j] == df_chain[i]['Residue'][j+1] - 1) and (df_chain[i]['ESA'][j+1]> 50):
                    
                        if df_chain[i]['pKa'][j] > df_chain[i]['pKa'][j+1]:
                            inter = 'non-reacting'
                        
                        else:
                        
                            if df_chain[i]['Amino Group'][j] == 'N':
                                inter = 'fast-reacting'
                            
                            else:
                            
                                if df_chain[i]['Secondary Structure'][j] == 'Helix':
                                    inter = 'slow-reacting'
                                
                                elif df_chain[i]['Secondary Structure'][j] == 'Strand':
                                    if df_chain[i]['pKa'][j] <= 10.3:
                                        if df_chain[i]['Area Of Lower Charge'][j] == 'Yes':
                                            inter = 'fast-reacting'
                                        else:
                                            inter = 'slow-reacting'
                                    
                                    else:
                                        inter = 'slow-reacting'
                                
                                else:
                                    if df_chain[i]['ESA'][j] >= 100:
                                        if df_chain[i]['Area Of Lower Charge'][j] == 'Yes':
                                            inter = 'fast-reacting'
                                    
                                        elif df_chain[i]['H-Donor'][j] == 'Yes':
                                            inter = 'fast-reacting'
                                    
                                        else:
                                            inter = 'slow-reacting'
                                        
                                    else:
                                        inter = 'slow-reacting'
                                        
                
                
                    else:
                    
                        if df_chain[i]['Amino Group'][j] == 'N':
                            inter = 'fast-reacting'
                            
                        else:
                            
                            if df_chain[i]['Secondary Structure'][j] == 'Helix':
                                inter = 'slow-reacting'
                                
                            elif df_chain[i]['Secondary Structure'][j] == 'Strand':
                                if df_chain[i]['pKa'][j] <= 10.3:
                                    if df_chain[i]['Area Of Lower Charge'][j] == 'Yes':
                                        inter = 'fast-reacting'
                                    else:
                                        inter = 'slow-reacting'
                                    
                                else:
                                    inter = 'slow-reacting'
                                
                            else:  # Coil, Turn, or other flexible structures
                                if df_chain[i]['ESA'][j] >= 100:
                                    if df_chain[i]['Area Of Lower Charge'][j] == 'Yes':
                                        inter = 'fast-reacting'
                                    
                                    elif df_chain[i]['H-Donor'][j] == 'Yes':
                                        inter = 'fast-reacting'
                                    
                                    else:
                                        inter = 'slow-reacting'
                                    
                                else:
                                    inter = 'slow-reacting'
                    
                                        
            interaction.append(inter)
        
    df['Interaction'] = interaction

    df.to_excel('Results.xlsx')

    return 

