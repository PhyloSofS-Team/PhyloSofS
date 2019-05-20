#!/usr/bin/env python
#coding: utf-8
import os, sys, subprocess, re, time
start = time.time()

with open(sys.argv[1].split("/")[-1][:-4]+".txt", "w+") as out:
    subprocess.call(["perl",
                      "/home/labeeuw/Documents/softwares/JPred_API_client_v_1_5/jpredapi",
                      "submit",
                      "mode=msa",
                      "format=msf",
                      #"email=antoine.labeeuw@etu.univ-nantes.fr",
                      "file="+sys.argv[1],
                      "name="+sys.argv[1].split("/")[-1][:-4]], stdout=out)
with open(sys.argv[1].split("/")[-1][:-4]+".txt", "r") as fileIn:
    for line in fileIn:
        if re.search("jobid:", line):
            jobid = line.split(" ")[-1]
b = subprocess.call(["perl",
                     "/home/labeeuw/Documents/softwares/JPred_API_client_v_1_5/jpredapi",
                     "status",
                     "jobid="+jobid,
                     "getResults=yes",
                     "checkEvery=5",
                     "silent"])
end = time.time()
print("time = ", end-start)
