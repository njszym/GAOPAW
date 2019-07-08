import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def read_values(col_index,dataset):
    col = df[:,col_index][4:]
    diff_list = []
    if dataset == 'GBRV':
        index = 1
    while index < len(col):
        diff_list.append(float(col[index]))
        index += 3
    return diff_list

elem_list = ['H','Li','Be','B','C','N','O','F','Na','Mg','Al','Si','P','S','Cl','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br']
elem_index_list = [i for i in range(len(elem_list))]

df = np.array(pd.read_csv('Accuracy.csv'))

lat_diff_list = read_values(1,'GBRV')
delta_diff_list = read_values(2,'GBRV')
gap_diff_list = read_values(3,'GBRV')
phon_diff_list = read_values(4,'GBRV')

lat_elem_list = []
lat_elem_index_list = []
new_lat_diff_list = []
for (elem,elem_index,lat_diff) in zip(elem_list,elem_index_list,lat_diff_list):
    if str(lat_diff) != 'nan':
        new_lat_diff_list.append(lat_diff)
        lat_elem_list.append(elem)
        lat_elem_index_list.append(elem_index)
lat_diff_list = new_lat_diff_list

delta_elem_list = []
delta_elem_index_list = []
new_delta_diff_list = []
for (elem,elem_index,delta_diff) in zip(elem_list,elem_index_list,delta_diff_list):
    if str(delta_diff) != 'nan':
        new_delta_diff_list.append(delta_diff)
        delta_elem_list.append(elem)
        delta_elem_index_list.append(elem_index)
delta_diff_list = new_delta_diff_list

gap_elem_list = []
gap_elem_index_list = []
new_gap_diff_list = []
for (elem,elem_index,gap_diff) in zip(elem_list,elem_index_list,gap_diff_list):
    if str(gap_diff) != 'nan':
        new_gap_diff_list.append(gap_diff)
        gap_elem_list.append(elem)
        gap_elem_index_list.append(elem_index)
gap_diff_list = new_gap_diff_list

phon_elem_list = []
phon_elem_index_list = []
new_phon_diff_list = []
for (elem,elem_index,phon_diff) in zip(elem_list,elem_index_list,phon_diff_list):
    if str(phon_diff) != 'nan':
        new_phon_diff_list.append(phon_diff)
        phon_elem_list.append(elem)
        phon_elem_index_list.append(elem_index)
phon_diff_list = new_phon_diff_list

fig, ax = plt.subplots(nrows=4,ncols=1)
fig.set_figwidth(14)
fig.set_figheight(18)

for i in range(4):
    ax[i].set_xticks(elem_index_list)
    ax[i].set_xticklabels(elem_list)
    ax[i].set_xlim(min(elem_index_list),max(elem_index_list))
    ax[i].set_ylabel('Average Absolute Error (%)',fontsize=16,labelpad=12.0)
    ax[i].tick_params(axis='x',labelsize=13.5)
    ax[i].tick_params(axis='y',labelsize=13)
    ax[i].minorticks_on()
    ax[i].grid(zorder=0)

ax[0].bar(lat_elem_index_list,lat_diff_list,zorder=3)
ax[0].set_ylim(0,max(lat_diff_list)+1)
ax[0].set_title('Lattice Constant',fontsize=20)

ax[1].bar(delta_elem_index_list,delta_diff_list,zorder=3)
ax[1].set_ylim(0,max(delta_diff_list)+1)
ax[1].set_title('Delta-Factor',fontsize=20)
ax[1].set_ylabel('Absolute Error',fontsize=16,labelpad=12.0)
ax[1].set_ylim(0,max(delta_diff_list)+0.5)

ax[2].bar(gap_elem_index_list,gap_diff_list,zorder=3)
ax[2].set_ylim(0,max(gap_diff_list)+3)
ax[2].set_title('Band Gap',fontsize=20)

ax[3].bar(phon_elem_index_list,phon_diff_list,zorder=3)
ax[3].set_ylim(0,max(phon_diff_list)+5)
ax[3].set_title('Phonon Frequency',fontsize=20)

plt.tight_layout()
plt.savefig('Accuracy.png',dpi=500)
plt.close()
