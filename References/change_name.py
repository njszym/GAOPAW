import os

for file in os.listdir('.'):
    if file[:2] == 'AE':
        os.rename(file,file[3:])
