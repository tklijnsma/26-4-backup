import itertools as it


jets = [ 'j'+str(i) for i in range(1,9) ]

qs = [ 'q'+str(i) for i in range(1,4) ]

links = [ [ 1, 5, 7 ], [ 1, 2, 6, 7 ], [ 4, 5 ] ]




print jets
print qs

print links



i_list = [ 0, 0, 0 ]

max_i_list = [ len(i)-1 for i in links ]


# Create unique lists of permutations

s0 = [ 0, 0, 1 ]

s_start = [ 0, 0, 0 ]
s = s_start


for count in range(3):

    for x in set(it.permutations(s)):

        for y in set(it.permutations(s0)):

            print [ i+j for (i,j) in zip(x,y) ]

    s0[0] += 1
        







"""
within_indices = True

for sumi in range( max( max_i_list ) ):



    for (i, i_max) in zip(i_list, max_i_list):
        if i > i_max:
            within_indices = False
"""
