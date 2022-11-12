# Import all COD files and generate all supercell files

import urllib.request


file = open('COD_ID_List.txt')
count = 0
for line in file:
    url = 'http://www.crystallography.net/cod/' + line.strip() + '.cif'
    filename = './CIF/' + line.strip() + '.cif'
    urllib.request.urlretrieve(url, filename)
    count += 1

print("Download " +str(count) + " CIF files success")