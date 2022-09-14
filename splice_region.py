#!/usr/bin/env python

def splice_region(region, introns):
    '''
    Take a sequence string and a list of lists containing
    intron borders and splice out the introns from the 
    sequence string
    '''
    # -1 from intron start (='exon' end)
    # +1 from intron end   (='exon' start)
    introns = [ [b[0]-1 , b[1]+1] for b in introns ]
        
    # flatten list and convert to 0-indexing
    f = [ x-1 for b in introns for x in b ]
    # add start and end (0-indexed)
    f.insert(0, 0)
    f.append(len(region)-1)

    # repack list
    borders = [ [f[i],f[i+1]] for i in range(0, len(f), 2) ]
    # construct spliced region
    spliced_region = ''
    for b in borders:
        spliced_region += region[ b[0]:b[1]+1 ]
    return spliced_region


seq = '1234567+++++891+++++23456'
coords = [ [8,12] , [16,20] ]

print( splice_region(seq, coords) )
