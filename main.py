import pandas as pd
import itertools

from SGC import *


Freqdf=pd.read_excel('AAfreqs.xlsx')

print(Freqdf)

Freqdf=Freqdf.set_index('amino')

Freqlookup=Freqdf.to_dict()

print(Freqlookup)

FreqVal=Freqlookup['p(a)'] #use FreqVal to get the frequecy of aa within p(a) key
print(FreqVal['Trp']) # returns 0.0109 correct value

Syndf=pd.read_excel('SynCodonsSGC.xlsx')
print(Syndf)
Syndf=Syndf.set_index('amino')
Synlookup=Syndf.to_dict()
print(Synlookup)

SynVal=Synlookup['n(a(c))'] # gets syn codons number from within n(a(c)) key
print(SynVal['Glu']) # returns 2 correct value



#%%
'''
for i in sgc.keys():
    print(SynVal[sgc[i]])
# comes back as sgc is not deifned, doing something obvious wrong?
'''


#%%


def get_gilis_val(df, column_label, row_label):
    '''
    This function returns a Gilis value
    :param df: dataframe
    :param column_label: amino code (string)
    :param row_label: amino code (string)
    :return: value (int)
    '''
    return df[column_label][row_label]

def compare_codon(a, b):
    '''
    Compare two codons and return True if they only differ by one base
    :param a: codon (string) e.g. 'UUU'
    :param b: codon (string) e.g. 'UUU'
    :return: boolean
    '''
    if len([i for i in range(len(a)) if a[i] != b[i]]) == 1:
        return True
    else:
        return False

#%%
        '''
SynCodonData=pd.read_excel("SynCodonsSGC.xlsx")

def get_SynCodons(df, row_label, column_label):
    df=SynCodonData
    return df[row_label][column_label]

get_SynCodons(df,sg)
'''
# so, itertools doesn't compare the same codons twice
# if you run the file you'll see the first valid comparison is UUU and UUC
# but, once UUU is exhausted and itertools moves to UUC, UUC is not then compared to UUU
# in theory this is fine for overall cot error, as all end values will be half that of-
#-what it would be if the same codons were compared again.
# this is becuase (currently) the gilis value would be the same for the secound comparison
# But, how does this affect N?
# since N is a normalisation for any c, 
# is N correct if the duplicate comparisons aren't made?
# if so, then N is even larger ( twice as large) then it is currently
# (this is notcounting anything else that needs to happen to N of course)
# so should I ue itertools.combinations_with_replacements ?
# actually i beive it will need to be, as once the amino acid frequency-
#- and the amount of synomounous codon d c is included the overall error cost-
# will be different, e.g. between UUA and UUU, currently only compared once-
# however when compared twice in the complete equation the UUA and UUU-
#- will have different frequency and synomous codn values since UUA codes Leu-
#- where as UUU encodes Phe
    



#%%

purines=['A','G']
pyrimidines=['C','U']

def compare_list_elements(df, sgc):
    sum = 0
    count = 0
    p=0
    for a, b in itertools.permutations(sgc.keys(), 2):
        if ('Ter' not in [sgc[a]]) and (a != b):
            count += 1
            if compare_codon(a, b):
                if a[0:1] != b[0:1]:
                    if (a[0:1] in purines) and (b[0:1] in purines):
                        p=(1/5.7)
                        aafreq= FreqVal[sgc[a]]
                        synnum= SynVal[sgc[a]]
                        gilis_value = get_gilis_val(df, sgc[a], sgc[b])
                        sum +=(aafreq/synnum)*(p*gilis_value)
                        print(count, [a, b], [sgc[a], sgc[b]], aafreq, synnum, p, gilis_value, sum)
                    elif (a[0:1] in pyrimidines) and (b[0:1] in pyrimidines):
                        p=(1/5.7)
                        aafreq= FreqVal[sgc[a]]
                        synnum= SynVal[sgc[a]]
                        gilis_value = get_gilis_val(df, sgc[a], sgc[b])
                        sum +=(aafreq/synnum)*(p*gilis_value)
                        print(count, [a, b], [sgc[a], sgc[b]], aafreq, synnum, p, gilis_value, sum)
                    else:
                        p=(0.5/5.7)
                        aafreq= FreqVal[sgc[a]]
                        synnum= SynVal[sgc[a]]
                        gilis_value = get_gilis_val(df, sgc[a], sgc[b])
                        sum +=(aafreq/synnum)*(p*gilis_value)
                        print(count, [a, b], [sgc[a], sgc[b]], aafreq, synnum, p, gilis_value, sum)
                elif a[1:2] != b[1:2]:
                    if (a[1:2] in purines) and (b[1:2] in purines):
                        p=(0.5/5.7)
                        aafreq= FreqVal[sgc[a]]
                        synnum= SynVal[sgc[a]]
                        gilis_value = get_gilis_val(df, sgc[a], sgc[b])
                        sum +=(aafreq/synnum)*(p*gilis_value)
                        print(count, [a, b], [sgc[a], sgc[b]], aafreq, synnum, p, gilis_value, sum)
                    elif (a[1:2] in pyrimidines) and (b[1:2] in pyrimidines):
                        p=(0.5/5.7)
                        aafreq= FreqVal[sgc[a]]
                        synnum= SynVal[sgc[a]]
                        gilis_value = get_gilis_val(df, sgc[a], sgc[b])
                        sum +=(aafreq/synnum)*(p*gilis_value)
                        print(count, [a, b], [sgc[a], sgc[b]], aafreq, synnum, p, gilis_value, sum)
                    else:
                        p=(0.1/5.7)
                        aafreq= FreqVal[sgc[a]]
                        synnum= SynVal[sgc[a]]
                        gilis_value = get_gilis_val(df, sgc[a], sgc[b])
                        sum +=(aafreq/synnum)*(p*gilis_value)
                        print(count, [a, b], [sgc[a], sgc[b]], aafreq, synnum, p, gilis_value, sum)
                elif a[2:3] != b[2:3]:
                    p=(1/5.7)
                    aafreq= FreqVal[sgc[a]]
                    synnum= SynVal[sgc[a]]
                    gilis_value = get_gilis_val(df, sgc[a], sgc[b])
                    sum +=(aafreq/synnum)*(p*gilis_value)
                    print(count, [a, b], [sgc[a], sgc[b]], aafreq, synnum, p, gilis_value, sum)
    return sum



# prototype of error cost function with amino acid frequencies-
#- and synonomous codons incorporated see 13 lines down
'''
def compare_list_elements(df, sgc):
    sum = 0
    count = 0
    p=0
    for a, b in itertools.permutations(sgc.keys(), 2):
        if ('Ter' not in [sgc[a], sgc[b]]) and (a != b):
            count += 1
            if compare_codon(a, b):
                if a[0:1] != b[0:1]:
                    if (a[0:1] in purines) and (b[0:1] in purines):
                        p=(1/5.7)
                        aminofreq=FreqVal[sgc[a]]
                        synnum=SynVal[sgc[a]]
                        gilis_value = get_gilis_val(df, sgc[a], sgc[b])
                        sum += (aminofrq/synnum) * (p*gilis_value) # new sum line
                        print(count, [a, b], [sgc[a], sgc[b]], gilis_value, sum)
                    elif (a[0:1] in pyrimidines) and (b[0:1] in pyrimidines):
                        p=(1/5.7)
                        gilis_value = get_gilis_val(df, sgc[a], sgc[b])
                        sum +=(p*gilis_value)
                        print(count, [a, b], [sgc[a], sgc[b]], gilis_value, sum)
                    else:
                        p=(0.5/5.7)
                        gilis_value = get_gilis_val(df, sgc[a], sgc[b])
                        sum +=(p*gilis_value)
                        print(count, [a, b], [sgc[a], sgc[b]], gilis_value, sum)
                elif a[1:2] != b[1:2]:
                    if (a[1:2] in purines) and (b[1:2] in purines):
                        p=(0.5/5.7)
                        gilis_value = get_gilis_val(df, sgc[a], sgc[b])
                        sum +=(p*gilis_value)
                        print(count, [a, b], [sgc[a], sgc[b]], gilis_value, sum)
                    elif (a[1:2] in pyrimidines) and (b[1:2] in pyrimidines):
                        p=(0.5/5.7)
                        gilis_value = get_gilis_val(df, sgc[a], sgc[b])
                        sum +=(p*gilis_value)
                        print(count, [a, b], [sgc[a], sgc[b]], gilis_value, sum)
                    else:
                        p=(0.1/5.7)
                        gilis_value = get_gilis_val(df, sgc[a], sgc[b])
                        sum +=(p*gilis_value)
                        print(count, [a, b], [sgc[a], sgc[b]], gilis_value, sum)
                elif a[2:3] != b[2:3]:
                    p=(1/5.7)
                    gilis_value = get_gilis_val(df, sgc[a], sgc[b])
                    sum +=(p*gilis_value)
                    print(count, [a, b], [sgc[a], sgc[b]], gilis_value, sum)
    return sum
'''

'''  
def compare_list_elements(df, sgc):
    sum = 0
    count = 0
    for a, b in itertools.combinations(sgc.keys(), 2):
        if ('Ter' not in [sgc[a], sgc[b]]) and (a != b):
            count += 1
            if compare_codon(a, b):
                gilis_value = get_gilis_val(df, sgc[a], sgc[b])
                sum += gilis_value
                print(count, [a, b], [sgc[a], sgc[b]], gilis_value, sum)
    return sum
'''

'''
purines=['A','G']
pyrimidines=['C','U']
'''
# need to run whole file or N  won't work
# N=5.7
# N may simpley be 5.7, wait for confrim, but loop still useable for p definition.

'''
def acquireN(sgc):
    N=0
    for a,b in itertools.combinations(sgc.keys(), 2):
        if('Ter' not in [sgc[a],sgc[b]]) and (a != b):
             if a[0:1] == b[0:1]:
                    if a[1:2] == b[1:2]:
                        if a[2:3] != b[2:3]:
                            N += 1
                    elif a[2:3] == b[2:3]:
                        if (a[1:2] in purines) and (b[1:2] in purines):
                            N += 0.5 #purine to purine transition
                        if (a[1:2] in pyrimidines) and (b[1:2] in pyrimidines):
                            N += 0.5 # pyrimidine to pyrimidine transition
                        else:
                            N += 0.1 # one purine, one pyrimidine = transversion
             elif a[1:2] == b[1:2]:
                    if a[2:3] == b[2:3]:
                        if (a[0:1] in purines) and (b[0:1] in purines):
                            N += 1
                        if (a[0:1] in pyrimidines) and (b[0:1] in pyrimidines):
                            N += 1
                        else:
                            N += 0.5
    return N
    print(N)                   
 '''   

# is the loop correct? i think so but N seems very large.
 
'''
acquireN(SGC)
'''

# is this right? it would make N pretty big. or, since its sum p(c|c')= 1 for any c -
#- does it need divided by 64,61(without stop) or 9, because each codon would only give-
# a p(c|c') value that wasn't Zero 9 times (each single base mutation)

#%%
# order to run function
if __name__ == "__main__":
    gilis_df = pd.read_excel('gilis.xlsx')
    sum = compare_list_elements(gilis_df, SGC)
    print(sum)


# sum p(c|c')= N intially make a function that gets sum
# then rerun but * gilis values by p where p is (correct numerator)/N
#N is the sum of p(c|c') numerators ( e.g 1 for a third base change)
#therefore need to make a function that runs through codon comparisons 
# then sums the p(c|c') values 
#not nesessacary for function to returh gilis values 
#after executeing the N function 
#then run gilis compariosn mutiplying the gilis values by th nesacary value over N









