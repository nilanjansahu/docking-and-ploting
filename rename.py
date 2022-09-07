import os

for filename in os.listdir(os.getcwd()):
    os.rename(filename, filename.replace(',','_'))