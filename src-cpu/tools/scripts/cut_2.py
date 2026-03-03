#!/usr/bin/env python3
import json
## cut the two pair.json and thermo.json together
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-t','--subdirs', nargs='*',help="thermo json file")
parser.add_argument("-p", "--pair", help="pair json file", action='store')
parser.add_argument("-o", "--out", help="out file dir and name", action='store')
parser.add_argument("-n", "--name",nargs='*', help="name of the data", action='store')
## 需要先将两个thermo结果合起来再去做其他的
args = parser.parse_args()

data=[]
with open(args.pair, 'r') as f1:
    data_pair = json.load(f1)
for subdir in args.subdirs:
        with open(subdir,'r') as f:
            data+=[json.load(f)]

keys=iter(args.name)
values=iter(data)
# {"pair": data_pair, key:value ,"type":args.name for key, value in zip(keys, values)}
with open(args.out, 'w') as f:
    # json.dump({"pair": data_pair, key:value ,"type": args.name for key, value in zip(keys, values)}, f, indent=4, separators=(',', ':'))
    json.dump({"pair":data_pair,"data":{key: value for key, value in zip(keys, values)},"type":args.name}, f, indent = 4, separators = (',', ':'))

    # json.dump({"pair":data_pair,"result":data,"type":args.name}, f,indent=4, separators=(',', ':'))
    f.write('\n')
