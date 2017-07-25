#!/usr/bin/env python 

import sys,argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='checks a pedigree against akt mendel output.')
    parser.add_argument('-fam', metavar='fam', type=str, help='plink fam file',required=True)
    parser.add_argument('-mendel', metavar='mendel', type=str, help='akt mendel output',required=True)
    parser.add_argument('-error', metavar='error', type=float, help='maximum mendel inconsistency rate',default = 0.05)
    args = parser.parse_args()

    d = {}
    for linenum,line in enumerate(open(args.mendel)):
        if linenum>0:
            ped_id,child_id,dad_id,mum_id,dad_gt,mum_gt,child_rr,child_ra,child_aa,nerror,error_rate,het_rate = line.split()
            if dad_id==".":
                dad_id="0"
            if mum_id==".":
                mum_id="0"

            key=(child_id,dad_id,mum_id)
            if key not in d:
                d[key] = [0,0]

            d[key][0] += int(nerror)
            if dad_gt!="RR" or mum_gt!="RR":
                d[key][1] += int(child_rr)
            d[key][1] += int(child_ra)
            d[key][1] += int(child_aa)

    for line in open(args.fam):
        fid,sampleid,dad,mum,sex = line.split()[:5]
        if dad!="0" or mum!="0":
            key = (sampleid,dad,mum)
            error = float(d[key][0])/float(d[key][1])
            if error<args.error:
                print line.strip()
            else:
                sys.stderr.write("WARNING: trio %s had a high mendel inconsistency rate (%f)\n"%(line.strip(),error))
        else:
            print line.strip()
