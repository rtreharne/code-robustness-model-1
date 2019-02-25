import pandas as pd
from collections import OrderedDict
from itertools import islice
import copy
import random
def get_ordered_gc(fname):

    # read in all data from blastocrithidia spreadsheet (block format)
    df = pd.read_excel(fname, header=None)

    # split columns into pairs, each a new dataframe with columns "codon" and "amino"
    df_col_pairs = [df.iloc[:, i*2:i*2+2].rename(columns={i*2:"codon", i*2+1:"amino"}) for i in range(0, int(len(df.columns.values)/2))]

    return pd.concat(df_col_pairs).set_index("codon").to_dict(into=OrderedDict)['amino']

class CodeOptimise:

    def __init__(self, GC):
        self.GC = GC
        self.blocks = self.block_code()
        print(self.blocks)
        self.pairwise_swaps()

    def block_code(self, block_size=4):
        blocks = [OrderedDict(islice(self.GC.items(), i*block_size, i*block_size+block_size)) for i in range(0, int(len(self.GC)/block_size))]
        return blocks

    def pairwise_swaps(self):

        print(self.swap(0, 1))
        print(self.swap(0, 2))
        print(self.swap(0,3))




    def swap(self, index_1, index_2):
        od_1_vals = list(self.blocks[index_1].values())
        od_2_vals = list(self.blocks[index_2].values())

        new_blocks = copy.deepcopy(self.blocks)

        for i, codon in enumerate(self.blocks[index_1]):
            new_blocks[index_1][codon] = od_2_vals[i]
        for i, codon in enumerate(self.blocks[index_2]):
            new_blocks[index_2][codon] = od_1_vals[i]

        return new_blocks


if __name__ == "__main__":
    gc = get_ordered_gc("../blastocrithidia_code.xlsx")
    opt = CodeOptimise(gc)


'''
print(opt.swap(0,1))
print(opt.swap(0,2)) 
print(opt.swap(1,3)) #swap is now working so prev swaps not in next print
'''

SwapList= [] # to assemble list of swapped dicts

# print(opt.blocks(0)) reurns 'list item not callable'

'''
for i,j in opt.blocks.items():
    print(i,j) # returns " 'list' object has no attribute 'items'"
'''

'''
for i in opt:
    for j in opt:
        SwappedDict={}
        if (i != j) and ('Glu/Stp' not in [[i],[j]]):
          SwappedDict = opt.swap(i,j)
          SwapList.append(SwappedDict)
          
print(SwapList) # returns CodeOptimise object is not interable
'''
'''
for i,j in enumerate(opt.blocks):
    SwappedDict= [] # tried {} as well
    if (i != j) and ('Glu/Stp' not in [[i], [j]]):
        SwappedDict = opt.swap(i,j)
        SwapList.append(SwappedDict)
        
print(SwapList) #returns list indices must be int or slices, not collections.OrderedDict
'''
'''
for i in opt.blocks:
    print(i)
print('i printed') # returns normal ordered dict
'''

'''
for (i) in opt.blocks:
    for (j) in opt.blocks:
        if (i !=j):
            NewBlocks=opt.swap(i,j)
            SwapList.append(NewBlocks)
print(SwapList) # returns list indices must be int or slice, not collections.ordererDict        
'''

'''
for (int(i) in range(0,16)) in opt.blocks:
    for (int(j) in range(0,16)) in opt.blocks:
        if (i != j):
            NewBlocks=opt.swap(i,j)
            SwapList.append(NewBlocks)
                     
print(SwapList) # returns can't assign to comparison (refering to i)
'''
print(' rand starts here')

# tried useing import random
'''
i = random.randint(0,16)
j = random.randint(0,16)
# if i and j outside loop doesn't work
'''
# opt.blocks(3) returns 'list oblect is not callable

count= 0
'''
for i in opt.blocks:
    for j in opt.blocks:
        i = random.randint(0,15)
        j = random.randint(0,15)
        if (i != j):
            NewBlocks= opt.swap(i,j)
            SwapList.append(NewBlocks)
            count += 1
            
            
print(SwapList)
print(count) # this is definitely doing something, might be whats needed
# prints alot so hard to tell, but think its doing swaps, but in a random order
# if rand order not important then this should work
# random order shouldn't matter for finding best code of all, but willl be akward
# count number is different at end of each run, all around 233-246 
# 91 true swaps should be made, but there will be repeats from this method.
# if so need to remove duplicates from SwapList
# if i had to gess this methd can't be relied upon to get all codon arrange ments
'''  
index1=list(range(0,16))
index2=list(range(0,16))

for i in index1:
    for j in index2:
        if (i != j):
            NewBlocks= opt.swap(i,j)
            SwapList.append(NewBlocks)
            count += 1   

print(SwapList)
print(count)
# reliably results in 240 count, witch is 16 *15,
#16 is the nmber of blocks, 15 th number of swaps each block can do
# this method will also have repats so needs to have repeats removed from Swaplist
# need to update to prevent swapping Glu/Stp and intially Trp
'''
SwapListNoDups= set(SwapList)

Swaps=list(SwapListNoDups)

# this was supposed to remove duplicates but didn't work

SwapsCount=0
    
for i in Swaps:
    print(i)
    SwapsCount +=1

print('Swap Count equals', SwapsCount) 
'''     