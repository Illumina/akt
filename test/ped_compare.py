import sys

def readPed(fname):
    ret = []
    for line in open(fname):
        ret.append( line.split()[1:4] )
    return(ret)

def pedSummary(ped):
    ntrio = len([val for val in ped if val[1]!="0" and val[2]!="0"])
    nduo = len([val for val in ped if val[1]!="0" or val[2]!="0"])
    nduo -= ntrio
    print ntrio,"trios"
    print nduo,"duos"

#extracts trios. orders parents alphabetically. returns set of "id-pid-mid"
def getTrios(ped):
#    import pdb;pdb.set_trace()
    ids = set([val[0] for val in ped])
    return(set( ["%s-%s-%s"%(val[0],min(val[1],val[2]),max(val[1],val[2])) for val in ped if val[1] in ids and val[2] in ids] ))

#extracts duos. orders alphabetically since direction of relationship cannot be ascertained.
def getDuos(ped):
#    import pdb;pdb.set_trace()
    ids = set([val[0] for val in ped])
    mums = set( ["%s-%s"%tuple(sorted([val[0],val[1]])) for val in ped if val[1] in ids and val[2] not in ids])
    dads = set( ["%s-%s"%tuple(sorted([val[0],val[2]])) for val in ped if val[2] in ids and val[1] not in ids])
    return(mums.union(dads))
    

if len(sys.argv)!=3:
    print "Usage: ped_compare.py test.fam truth.fam"
else:
    test_ped=readPed(sys.argv[1])
    pedSummary(test_ped)
    truth_ped=readPed(sys.argv[2])
    pedSummary(truth_ped)
    ids1 = set( [val[0] for val in test_ped] )
    ids2 = set( [val[0] for val in truth_ped] )
    ids = ids1.intersection(ids2)
    print len(ids),"ids intersecting"

    test_ped = [val for val in test_ped if val[0] in ids]
    truth_ped = [val for val in truth_ped if val[0] in ids]
    trios1 = getTrios(test_ped)
    trios2 = getTrios(truth_ped)
    # print trios1
    # print trios2
    print "\nTRIOS:"
    print "Exact match?",trios1==trios2
    print "FP:",trios1.difference(trios2)
    print "FN:",trios2.difference(trios1)

    duos1 = getDuos(test_ped)
    duos2 = getDuos(truth_ped)
    print "\nDUOS:"
    print "Exact match?",duos1==duos2
    print "FP:",duos1.difference(duos2)
    print "FN:",duos2.difference(duos1)
