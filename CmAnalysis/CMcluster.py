#!/bin/python

import argparse
import numpy as np

#Command line Argument Parsing
p=argparse.ArgumentParser(description='Parse SOM output and sort into clusters')
p.add_argument('SOM', help='R SOM bins file')
p.add_argument('-o', type=str, help='Output file name', default='clusters')
args = p.parse_args()

bins = {}
gbid_loci = {}
species_mapping = {}

with open(args.SOM, 'r') as SOM_OUT:
    ids = SOM_OUT.readline().split()

    for line in SOM_OUT:
        tmp = line.split()
        cluster = tmp[0]

        if cluster not in bins:
            bins[cluster] = []
            bins[cluster].append(tmp[1:])
        else:
            bins[cluster].append(tmp[1:])

        som_bin = tmp[1]
        vir_name = tmp[2].split('|')[0]
        vir_gbid = tmp[2].split('|')[1]
        loc_pos = tmp[2].split('|')[2]

        if vir_gbid not in gbid_loci:
            gbid_loci[vir_gbid] = []
            gbid_loci[vir_gbid].append((loc_pos, cluster))
        else:
            gbid_loci[vir_gbid].append((loc_pos, cluster))

        if vir_name not in species_mapping:
            species_mapping[vir_name] = []
            species_mapping[vir_name].append(vir_gbid)
        else:
            if vir_gbid not in species_mapping[vir_name]:
                species_mapping[vir_name].append(vir_gbid)
            else:
                pass

#we want to know which Genome has how many hits
#also we want to know where those hits are located
#And what type of loci we are seeing
with open(args.o+'_comprehensive', 'w') as SUMM:
    for v_species in species_mapping:
        SUMM.write(v_species+'\n')
        for v_gbid in species_mapping[v_species]:
            SUMM.write('\t'+v_gbid+'\t')
            for locus in gbid_loci[v_gbid]:
                SUMM.write(locus[0]+' '+locus[1]+'\t')
            SUMM.write('\n')


#Write Our Cluster information
with open(args.o, 'w') as OUT:
    print ids
    OUT.write('\t'.join(ids)+'\n')

    for clu in bins:
        OUT.write(clu+'\n')
        for hit in bins[clu]:
            OUT.write('\t')
            OUT.write(hit[0])
            OUT.write('\t')
            OUT.write("{:100s}".format(hit[1]))
            OUT.write('\t')
            for i, val in enumerate(hit[2:], 2): #Start at pos 2, enumerate not needed
                OUT.write("{:5.3f}".format(float(val))+'\t')
            OUT.write('\n')


#Print Cluster summary vectors
with open(args.o+'_summary', 'w') as OUT:
    for clu in bins:
        print 'Cluster ID: ', clu
        clu_summary = [[] for x in range(len(bins[clu][0][2:]))] #Matrix for values
        clu_virnames = {}
        clu_vircount = {}
        virname_gbIDs = {} # virus_name -> [gbID, gbID, gbID ...] e.g. Dengue_1 -> ['AF9233', 'DF4345G3',...]

        for hit_vector in bins[clu]:
            virname, vir_id = hit_vector[1].split('|')[:2]

            #Count
            if virname not in clu_virnames:
                clu_virnames[virname] = 1
                virname_gbIDs[virname] = []
                virname_gbIDs[virname].append(vir_id)
            else:
                clu_virnames[virname] += 1
                virname_gbIDs[virname].append(vir_id)

            if vir_id not in clu_vircount:
                clu_vircount[vir_id] = 1
            else:
                clu_vircount[vir_id] += 1

            #Append E-value Data
            for i, val in enumerate(hit_vector[2:]):
                clu_summary[i].append(float(val))

        #Calc Mean/Stddev
        means = []
        stddevs = []
        for i, vals in enumerate(clu_summary):
            data = np.array(vals, dtype=float)
            means.append(data.mean())
            stddevs.append(data.std())

        #Gather information about hit multiplicity
        multiplicity = {}
        for virname in virname_gbIDs:
            hitcounter = {}
            for gb_id in virname_gbIDs[virname]:
                hits = clu_vircount[gb_id]  #Nr of loci in this genome
                if hits not in hitcounter:
                    hitcounter[hits] = 1
                else:
                    hitcounter[hits] += 1


            multiplicity[virname] = hitcounter

        #Output
        OUT.write(clu+'\n\t\t')                #Cluster ID
        OUT.write('\t'.join(ids[3:])+'\n')     #CM names
        OUT.write('\t\t')
        for m in means:
            OUT.write("{:5.3f}".format(float(m))+'\t')  #CM mean E-values
        OUT.write('\n')
        OUT.write('\t\t')
        for m in stddevs:
            OUT.write("{:5.3f}".format(float(m))+'\t')  #CM E-values Stddevs
        OUT.write('\n\n')
        for virname in clu_virnames:
            OUT.write('\t'+virname+'\t'+str(clu_virnames[virname])+'\t\t')
            for m in multiplicity[virname]:
                OUT.write(str(m)+' Hits: '+str(multiplicity[virname][m])+'\t')
            OUT.write('\n')
        OUT.write('\n')
