#!/usr/bin/env python3

'''
Parallel neoantigen relatived shell
'''

import os, sys
import pymongo
import re
import socket
import glob

pdir = sys.argv[1]
prj = sys.argv[2]
stype = sys.argv[3]
sample = sys.argv[4]

pepfile = []
pepnum = 0

pepfile=glob.glob(pdir+'/*.peptide')

pepnum=len(pepfile)

#for pfile in os.popen('ls ' + pdir + '/*.peptide').read().split("\n")[:-1]:
#    pepfile.append(pfile)
#    pepnum += 1

IP = socket.gethostbyname(socket.gethostname())
if re.match("192\.168\.1\.", IP):
    IP = "192.168.1.211"
else:
    IP = "192.168.11.101"
con = pymongo.MongoClient(IP, 27017)
db = con.pipeline
prjobj = db[prj]

'''
if pepnum == 0:
    prjobj.remove({"Name":"capneo_anno_" + sample})
    prjobj.remove({"Name":"netMHCIIpan_anno_" + sample})
'''

filearr = " ".join(pepfile)
length = ['8', '9', '10', '11']
for l in length:
    net1 = 'netMHCpan_' + stype + '_' + l + '_' + sample
    tap = 'capneo_tap_' + stype + '_' + l + '_' + sample
    if pepnum == 0:
        '''
        if stype == 'MT' and l == '8':
            prjobj.update({"Name":net1}, {"$set":{"Program":"template/mhc/netMHC.none.sh", "Output":"mhc/MHCIanno/" + sample, "Info":"final.report.txt"}})
        else:
            prjobj.remove({"Name":net1})
        prjobj.remove({"Name":tap})
        '''
        prjobj.update({"Name":net1}, {"$set":{"Program":"template/mhc/netMHCpan.none.sh", "Info":""}})
        prjobj.update({"Name":net1, "ArrayJob":{"$exists":True}}, {"$unset":{"ArrayJob":""}})
        prjobj.update({"Name":tap}, {"$set":{"Program":"template/mhc/TAP.none.sh", "Info":""}})
        prjobj.update({"Name":tap, "ArrayJob":{"$exists":True}}, {"$unset":{"ArrayJob":""}})
    elif pepnum <= 5:
        prjobj.update({"Name":net1}, {"$set":{"Program":"template/mhc/netMHCpan.small.sh", "Info":filearr, "ArrayJob":pepnum}})
        prjobj.update({"Name":tap}, {"$set":{"Program":"template/mhc/TAP.small.sh", "Info":filearr, "ArrayJob":pepnum}})
    elif pepnum <= 15:
        prjobj.update({"Name":net1}, {"$set":{"Program":"template/mhc/netMHCpan.medium.sh", "Info":filearr, "ArrayJob":pepnum}})
        prjobj.update({"Name":tap}, {"$set":{"Program":"template/mhc/TAP.medium.sh", "Info":filearr, "ArrayJob":pepnum}})
    elif pepnum <= 50:
        prjobj.update({"Name":net1}, {"$set":{"Program":"template/mhc/netMHCpan.large.sh", "Info":filearr, "ArrayJob":pepnum}})
        prjobj.update({"Name":tap}, {"$set":{"Program":"template/mhc/TAP.large.sh", "Info":filearr, "ArrayJob":pepnum}})
    else:
        prjobj.update({"Name":net1}, {"$set":{"Program":"template/mhc/netMHCpan.super.sh", "Info":filearr, "ArrayJob":pepnum}})
        prjobj.update({"Name":tap}, {"$set":{"Program":"template/mhc/TAP.super.sh", "Info":filearr, "ArrayJob":pepnum}})
length = ['9', '10', '11', '12', '13', '14', '15']
for l in length:
    net2 = 'netMHCIIpan_' + stype + '_' + l + '_' + sample
    if pepnum == 0:
        '''
        if stype == 'MT' and l == '9':
            prjobj.update({"Name":net2}, {"$set":{"Program":"template/mhc/netMHC.none.sh", "Output":"mhc/MHCIIanno/" + sample, "Info":"transcript.fp.hotspot.dnavaf.CCLEexpStatus"}})
        else:
            prjobj.remove({"Name":net2})
        '''
        prjobj.update({"Name":net2}, {"$set":{"Program":"template/mhc/netMHCIIpan.none.sh", "Info":""}})
        prjobj.update({"Name":net2, "ArrayJob":{"$exists":True}}, {"$unset":{"ArrayJob":""}})
    elif pepnum <= 5:
        prjobj.update({"Name":net2}, {"$set":{"Program":"template/mhc/netMHCIIpan.small.sh", "Info":filearr, "ArrayJob":pepnum}})
    elif pepnum <= 15:
        prjobj.update({"Name":net2}, {"$set":{"Program":"template/mhc/netMHCIIpan.medium.sh", "Info":filearr, "ArrayJob":pepnum}})
    elif pepnum <= 50:
        prjobj.update({"Name":net2}, {"$set":{"Program":"template/mhc/netMHCIIpan.large.sh", "Info":filearr, "ArrayJob":pepnum}})
    else:
        prjobj.update({"Name":net2}, {"$set":{"Program":"template/mhc/netMHCIIpan.super.sh", "Info":filearr, "ArrayJob":pepnum}})
