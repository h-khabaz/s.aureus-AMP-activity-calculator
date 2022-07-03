#!/usr/bin/python2

__author__ = "Lukasz Pawel Kozlowski"
__email__  = "lukaszkozlowski.lpk@gmail.com"
__copyrights__ = "Lukasz Pawel Kozlowski"
__webserver__ = "http://isoelectric.org"
__license__ = "http://isoelectric.org/license.txt"

import sys
import os

# to ad new pKa sets just add new dictionary       
scales = {
"EMBOSS":     {'Cterm': 3.6, 'pKAsp': 3.9,  'pKGlu': 4.1, 'pKCys': 8.5, 'pKTyr': 10.1, 'pk_his': 6.5, 'Nterm': 8.6, 'pKLys': 10.8, 'pKArg': 12.5},
"DTASelect":  {'Cterm': 3.1, 'pKAsp': 4.4,  'pKGlu': 4.4, 'pKCys': 8.5, 'pKTyr': 10.0, 'pk_his': 6.5, 'Nterm': 8.0, 'pKLys': 10.0, 'pKArg': 12.0},
"Solomon":    {'Cterm': 2.4, 'pKAsp': 3.9,  'pKGlu': 4.3, 'pKCys': 8.3, 'pKTyr': 10.1, 'pk_his': 6.0, 'Nterm': 9.6, 'pKLys': 10.5, 'pKArg': 12.5}, 
"Sillero":    {'Cterm': 3.2, 'pKAsp': 4.0,  'pKGlu': 4.5, 'pKCys': 9.0, 'pKTyr': 10.0, 'pk_his': 6.4, 'Nterm': 8.2, 'pKLys': 10.4, 'pKArg': 12.0},
"Rodwell":    {'Cterm': 3.1, 'pKAsp': 3.68, 'pKGlu': 4.25,'pKCys': 8.33,'pKTyr': 10.07,'pk_his': 6.0, 'Nterm': 8.0, 'pKLys': 11.5, 'pKArg': 11.5},
"Patrickios": {'Cterm': 4.2, 'pKAsp': 4.2,  'pKGlu': 4.2, 'pKCys': 0.0, 'pKTyr':  0.0, 'pk_his': 0.0, 'Nterm': 11.2,'pKLys': 11.2, 'pKArg': 11.2},
"Wikipedia":  {'Cterm': 3.65,'pKAsp': 3.9,  'pKGlu': 4.07,'pKCys': 8.18,'pKTyr': 10.46,'pk_his': 6.04,'Nterm': 8.2, 'pKLys': 10.54,'pKArg': 12.48},
"Grimsley":   {'Cterm': 3.3, 'pKAsp': 3.5,  'pKGlu': 4.2, 'pKCys': 6.8, 'pKTyr': 10.3, 'pk_his': 6.6, 'Nterm': 7.7, 'pKLys': 10.5, 'pKArg': 12.04},
'Lehninger':  {'Cterm': 2.34,'pKAsp': 3.86, 'pKGlu': 4.25,'pKCys': 8.33,'pKTyr': 10.0, 'pk_his': 6.0, 'Nterm': 9.69,'pKLys': 10.5, 'pKArg': 12.4},
'Bjellqvist': {'Cterm': 3.55,'pKAsp': 4.05, 'pKGlu': 4.45,'pKCys': 9.0, 'pKTyr': 10.0, 'pk_his': 5.98,'Nterm': 7.5, 'pKLys': 10.0, 'pKArg': 12.0},   
'IPC_peptide':{'Cterm': 2.383, 'pKAsp': 3.887, 'pKGlu': 4.317, 'pKCys': 8.297, 'pKTyr': 10.071, 'pk_his': 6.018, 'Nterm': 9.564, 'pKLys': 10.517, 'pKArg': 12.503},    # IPC peptide
'IPC_protein':{'Cterm': 2.869, 'pKAsp': 3.872, 'pKGlu': 4.412, 'pKCys': 7.555, 'pKTyr': 10.85,  'pk_his': 5.637, 'Nterm': 9.094, 'pKLys': 9.052,  'pKArg': 11.84},     # IPC protein 
'Toseland':   {'Cterm': 3.19,'pKAsp': 3.6,  'pKGlu': 4.29,'pKCys': 6.87,'pKTyr': 9.61, 'pk_his': 6.33,'Nterm': 8.71, 'pKLys': 10.45, 'pKArg':  12},
'Thurlkill':  {'Cterm': 3.67,'pKAsp': 3.67, 'pKGlu': 4.25,'pKCys': 8.55,'pKTyr': 9.84, 'pk_his': 6.54,'Nterm': 8.0, 'pKLys': 10.4, 'pKArg': 12.0},
'Nozaki':     {'Cterm': 3.8, 'pKAsp': 4.0,  'pKGlu': 4.4, 'pKCys': 9.5, 'pKTyr': 9.6,  'pk_his': 6.3, 'Nterm': 7.5, 'pKLys': 10.4, 'pKArg': 12},   
'Dawson':     {'Cterm': 3.2, 'pKAsp': 3.9,  'pKGlu': 4.3, 'pKCys': 8.3, 'pKTyr': 10.1, 'pk_his': 6.0, 'Nterm': 8.2, 'pKLys': 10.5, 'pKArg':  12},   
          }

aaDict = {'Asp':'D', 'Glu':'E', 'Cys':'C', 'Tyr':'Y', 'His':'H', 
          'Lys':'K', 'Arg':'R', 'Met':'M', 'Phe':'F', 'Leu':'L', 
          'Val':'V', 'Ala':'A', 'Gly':'G', 'Gln':'Q', 'Asn':'N',
          'Ile':'I', 'Trp':'W', 'Ser':'S', 'Thr':'T', 'Sec':'U',
          'Pro':'P', 'Xaa':'X', 'Sec':'U', 'Pyl':'O', 'Asx':'B',
          'Xle':'J', }

acidic = ['D', 'E', 'C', 'Y']
basic = ['K', 'R', 'H']

pKcterminal = {'D': 4.55, 'E': 4.75} 
pKnterminal = {'A': 7.59, 'M': 7.0, 'S': 6.93, 'P': 8.36, 'T': 6.82, 'V': 7.44, 'E': 7.7} 

#N-teminus, middle, C-terminus
promost ={
'K':[10.00,  9.80, 10.30],
'R':[11.50, 12.50, 11.50],
'H':[ 4.89,  6.08,  6.89],
'D':[ 3.57,  4.07,  4.57],
'E':[ 4.15,  4.45,  4.75],
'C':[ 8.00,  8.28,  9.00],
'Y':[ 9.34,  9.84, 10.34],
'U':[ 5.20,  5.43,  5.60], # ref (http://onlinelibrary.wiley.com/doi/10.1002/bip.21581/pdf)
}
           
promost_mid = {
"G":[7.50, 3.70],
"A":[7.58, 3.75],
"S":[6.86, 3.61],
"P":[8.36, 3.40],
"V":[7.44, 3.69],
"T":[7.02, 3.57],
"C":[8.12, 3.10],
"I":[7.48, 3.72],
"L":[7.46, 3.73],
"J":[7.46, 3.73],
"N":[7.22, 3.64],
"D":[7.70, 3.50],
"Q":[6.73, 3.57],
"K":[6.67, 3.40],
"E":[7.19, 3.50],
"M":[6.98, 3.68],
"H":[7.18, 3.17],
"F":[6.96, 3.98],
"R":[6.76, 3.41],
"Y":[6.83, 3.60],
"W":[7.11, 3.78],
"X":[7.26, 3.57],   #avg
"Z":[6.96, 3.535],  #("E"+"Q")/2
'B':[7.46, 3.57],   #("N"+"D")/2
'U':[5.20, 5.60], 
'O':[7.00, 3.50],     
}

sample_protein_sequence = 'MKKMQSIVLALSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDHGWWKQHYEWRGNRWHLHGPPPPPRHHKKAPHDHHGGHGPGKHHR'

def predict_isoelectric_point_ProMoST(seq):
    '''Calculate isoelectric point using ProMoST model'''
    NQ = 0.0
    pH = 6.51             #starting po pI = 6.5 - theoretically it should be 7, but average protein pI is 6.5 so we increase the probability of finding the solution
    pHprev = 0.0         
    pHnext = 14.0        
    E = 0.01             #epsilon means precision [pI = pH +- E]
    temp = 0.01
    while 1:
            if seq[0] in promost.keys():   QN1=-1.0/(1.0+pow(10,(promost[seq[0]][2]-pH)))
            else: QN1=-1.0/(1.0+pow(10,(promost_mid[seq[0]][1]-pH)))
            #print 
            if seq[-1] in promost.keys():  QP2= 1.0/(1.0+pow(10,(pH-promost[seq[-1]][0])))
            else: QP2=1.0/(1.0+pow(10,(pH-promost_mid[seq[-1]][0])))
            
            QN2=-seq.count('D')/(1.0+pow(10,(promost['D'][1]-pH)))           
            QN3=-seq.count('E')/(1.0+pow(10,(promost['E'][1]-pH)))           
            QN4=-seq.count('C')/(1.0+pow(10,(promost['C'][1]-pH)))           
            QN5=-seq.count('Y')/(1.0+pow(10,(promost['Y'][1]-pH)))        
            QP1= seq.count('H')/(1.0+pow(10,(pH-promost['H'][1])))                            
            QP3= seq.count('K')/(1.0+pow(10,(pH-promost['K'][1])))           
            QP4= seq.count('R')/(1.0+pow(10,(pH-promost['R'][1])))                
        
            NQ=QN1+QN2+QN3+QN4+QN5+QP1+QP2+QP3+QP4  
#%%%%%%%%%%%%%%%%%%%%%%%%%   BISECTION   %%%%%%%%%%%%%%%%%%%%%%%%
            if NQ<0.0:              #we are out of range, thus the new pH value must be smaller                     
                    temp = pH
                    pH = pH-((pH-pHprev)/2.0)
                    pHnext = temp
                    #print "pH: ", pH, ", \tpHnext: ",pHnext
            else:
                    temp = pH
                    pH = pH + ((pHnext-pH)/2.0)
                    pHprev = temp
                    #print "pH: ", pH, ",\tpHprev: ", pHprev

            if (pH-pHprev<E) and (pHnext-pH<E): #terminal condition, finding pI with given precision
                    return pH               
        

def predict_isoelectric_point(seq, scale='IPC_protein'):
    """calculate isoelectric point using 9 pKa set model"""
    pKCterm = scales[scale]['Cterm']
    pKAsp = scales[scale]['pKAsp']
    pKGlu = scales[scale]['pKGlu']
    pKCys = scales[scale]['pKCys']
    pKTyr = scales[scale]['pKTyr']
    pKHis = scales[scale]['pk_his']
    pKNterm = scales[scale]['Nterm']
    pKLys = scales[scale]['pKLys'] 
    pKArg = scales[scale]['pKArg']
    pH = 6.51             #starting po pI = 6.5 - theoretically it should be 7, but average protein pI is 6.5 so we increase the probability of finding the solution
    pHprev = 0.0         
    pHnext = 14.0        
    E = 0.01             #epsilon means precision [pI = pH +- E]
    temp = 0.01
    nterm=seq[0]
    if scale=='Bjellqvist':
        if nterm in pKnterminal.keys():
            pKNterm = pKnterminal[nterm]
    
    cterm=seq[-1]
    if scale=='Bjellqvist':
        if cterm in pKcterminal.keys():
            pKCterm = pKcterminal[cterm] 
            
    while 1:             #the infinite loop
        QN1=-1.0/(1.0+pow(10,(pKCterm-pH)))                                        
        QN2=-seq.count('D')/(1.0+pow(10,(pKAsp-pH)))           
        QN3=-seq.count('E')/(1.0+pow(10,(pKGlu-pH)))           
        QN4=-seq.count('C')/(1.0+pow(10,(pKCys-pH)))           
        QN5=-seq.count('Y')/(1.0+pow(10,(pKTyr-pH)))        
        QP1=seq.count('H')/(1.0+pow(10,(pH-pKHis)))            
        QP2=1.0/(1.0+pow(10,(pH-pKNterm)))                
        QP3=seq.count('K')/(1.0+pow(10,(pH-pKLys)))           
        QP4=seq.count('R')/(1.0+pow(10,(pH-pKArg)))            
        NQ=QN1+QN2+QN3+QN4+QN5+QP1+QP2+QP3+QP4
        #print NQ
        #%%%%%%%%%%%%%%%%%%%%%%%%%   BISECTION   %%%%%%%%%%%%%%%%%%%%%%%%
        if NQ<0.0:              #we are out of range, thus the new pH value must be smaller                     
            temp = pH
            pH = pH-((pH-pHprev)/2.0)
            pHnext = temp
            #print "pH: ", pH, ", \tpHnext: ",pHnext
        else:
            temp = pH
            pH = pH + ((pHnext-pH)/2.0)
            pHprev = temp
            #print "pH: ", pH, ",\tpHprev: ", pHprev

        if (pH-pHprev<E) and (pHnext-pH<E): #terminal condition, finding pI with given precision
            return pH
                    
                    
def ipc_author_information():
    '''add information about IPC'''
    print( '==============================================================================================\n' )
    print( '\t\t\t\tIPC - ISOELECTRIC POINT CALCULATOR\n' )
    print( 'AUTHOR: \tLukasz Pawel Kozlowski, lukaszkozlowski.lpk@gmail.com' )
    print( 'COPYRIGHTS: \tLukasz Pawel Kozlowski\n' )
    print( 'WEB SERVER: \thttp://isoelectric.org' )
    print( 'LICENSE: \thttp://isoelectric.org/license.txt\n' )
    print( '\t\t\t\t\tJanuary 2016' )
    print( '==============================================================================================\n' )    

def error_information():
    '''information how to run IPC script'''
    
    print( "Usage: ipc <fasta_file> <pKa set> <output_file> <plot_file>\n" )
    info_string = '''<fasta_file>    protein sequence(s) in fasta format, see ./examples
<pKa set>       one from pKa sets which will be used to calculate pI, default 'ALL' 
                (report pI using all models), valid options are:
                'ALL', 'IPC_protein', 'IPC_peptide', 'Bjellqvist', 'Dawson', 'Grimsley', 
                'Toseland', 'EMBOSS', 'Kozlowski', 'DTASelect', 'Wikipedia', 'Rodwell', 
                'Patrickios', 'Sillero', 'Thurlkill', 'Solomon', 'Nozaki', 
                'Lehninger', 'ProMoST'
                
<output_file>   output of the program with pI predicted using selected model(s), default name 
                <fasta_file>.pI.txt
<plot_file>     virtual 2D-PAGE scatter plot (molecular weight vs. isoelectric point) 
                represented as heat map, this option is available only if numpy and matplotlib 
                and scipy are installed'''             
    print(info_string)    
        
    print( "\nE.g. python ipc.py NC_010473_Ecoli.faa" )
    print( "ipc ./examples/NC_010473_Ecoli.faa ALL out.txt out.png (if installed by setup.py)" )

    sys.exit(1)

def check_additional_libraries():
    ''' check libraries for plotting '''
    try: 
        import numpy as np
    except:
        print('Warning: numpy is missing, you will not get nice scatter plot, but program will produce the isoelectric point predictions anyway')
        print('         To install numpy on Ubuntu, just write in the terminal: sudo apt-get install python-numpy')
        return 0
    try: 
        import matplotlib.pyplot as plt
    except:
        print('Warning: matplotlib is missing, you will not get nice scatter plot, but program will produce the isoelectric point predictions anyway')
        print('         to install matplotlib on Ubuntu, just write in the terminal: sudo apt-get install python-matplotlib')
        return 0
    try: 
        from scipy.stats import gaussian_kde
    except:
        print('Warning: scipy is missing, you will not get nice scatter plot, but program will produce the isoelectric point predictions anyway')
        print('         to install scipy on Ubuntu, just write in the terminal: sudo apt-get install python-scipy')    
        return 0
    return 1    
        
def fasta_reader(fasta_string):
    """
    reads fasta file and return table [ [head1, seq1], [head2, seq2], ...]
    it is endure for most of errors like: multiple line for sequence, white spaces etc.
    """
    fasta_tab = ['>'+n+'\n' for n in ('\n'+fasta_string).split('\n>')[1:] ]
    fasta_list = []
    seq_tab = []
    for fas in fasta_tab:
        #solves problems with MAC, Linux and Windows new line characters
        fas = fas.replace('\r', '') 
        tmp = fas.split('\n')
        head = tmp[0]
        #deals with multiple lines in sequence
        seq = ''.join(tmp[1:]).upper().replace('*','')
        fasta_list.append([head, seq])
        seq_tab.append(seq)
    return fasta_list

def calculate_molecular_weight(seq):
    """molecular weight"""
    massDict = {'D':115.0886, 'E': 129.1155, 'C': 103.1388, 'Y':163.1760, 'H':137.1411, 
           'K':128.1741, 'R': 156.1875,  'M': 131.1926, 'F':147.1766, 'L':113.1594,
           'V':99.1326,  'A':  71.0788,  'G':  57.0519, 'Q':128.1307, 'N':114.1038, 
           'I':113.1594, 'W': 186.2132,  'S':  87.0782, 'P': 97.1167, 'T':101.1051, 
           'U':141.05,   'h2o':18.01524, 'X':        0, 'Z':128.6231, 'O':255.31, 
           'B':114.5962, 'J': 113.1594,}
    molecular_weight = massDict['h2o']
    for aa in seq:
            molecular_weight+=massDict[aa]       
    return molecular_weight
    
def make_heat_map(mw_tab, pI_tab, plot_file, input_pKa_set):
    """
    virtual 2D-PAGE scatter plot, heat map
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.stats import gaussian_kde
    if input_pKa_set=='ALL': input_pKa_set = 'IPC_protein'
    mw_tab = [n/1000 for n in mw_tab]
    x = pI_tab
    y = mw_tab

    # Calculate the point density
    xy = np.vstack([x,y])
    try: z = gaussian_kde(xy)(xy)
    except: return

    fig, ax = plt.subplots()
    ax.scatter(x, y, c=z, s=20, edgecolor='')
    plt.xlim(1, 14)
    plt.ylim(0, max(y))
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    ax.set_xlabel('Isoelectric point (%s)'%input_pKa_set, fontsize=22)
    ax.set_ylabel('Molecular weight (kDa)', fontsize=22)
    print('Plot file: '+plot_file)
    #plt.show()
    plt.savefig(plot_file)

if __name__ == '__main__': 
    ipc_author_information()
    try: 
        fasta_file = sys.argv[1].strip()
    except:
        error_information()
        
    available_pKa_sets = list(scales.keys())
    available_pKa_sets.append('ProMoST')
    available_pKa_sets.sort()
    
    try: 
        input_pKa_set = sys.argv[2].strip()
        if input_pKa_set not in available_pKa_sets:
            if input_pKa_set!='ALL':
                print("""Warning: provided pKa set "%s" is not valid, instead default "ALL" will be used
for more information run program without arguments: python ipc.py\n"""%input_pKa_set)
                input_pKa_set = 'ALL'
    except:
        input_pKa_set = 'ALL'
    
    try:
        output_file= sys.argv[3].strip()
    except:
        output_file = fasta_file+'.pI.txt'

    try:
        plot_file= sys.argv[4].strip()
    except:
        plot_file = fasta_file+'.png' 
    
    plot_libraries = check_additional_libraries()
    
    #read fasta file
    try: 
        fasta_tab = fasta_reader(open(fasta_file).read())
    except:
        print('Something wrong with fasta file')
        sys.exit(1)
        
    mw_tab = []
    pI_tab = []
    
    if input_pKa_set=='ALL': new_fasta_with_pI = '#'+', '.join(available_pKa_sets)+', Avg_pI'+os.linesep+os.linesep
    else: new_fasta_with_pI = '#'+input_pKa_set+os.linesep
        
    for query in fasta_tab:
        header, sequence = query
        #print(header+'\n'+sequence)
        
        mw = calculate_molecular_weight(sequence)
        mw_tab.append(mw)
        
        #extend header by molecular weight information
        header += '||Molecular weight: %s Da'%round(mw,2)
        
        if input_pKa_set!='ALL':
            if input_pKa_set == 'ProMoST':
                pI = predict_isoelectric_point_ProMoST(sequence)
            else:
                pI = predict_isoelectric_point(sequence, input_pKa_set)
            pI_tab.append(pI)
            new_fasta_with_pI += header+os.linesep+str(round(pI,3))+os.linesep+sequence+os.linesep
            
        else:
            #this is 'ALL' case, need to handle few things
            avg_ip_tab = []
            ip_tab_tmp = [] 
            for input_pKa_set_tmp in available_pKa_sets:
                #print(input_pKa_set)
                #print('\n')
                if input_pKa_set_tmp == 'ProMoST':
                    pI = predict_isoelectric_point_ProMoST(sequence)
                else:
                    pI = predict_isoelectric_point(sequence, input_pKa_set_tmp)
                
                #this one is needed for plot
                if input_pKa_set_tmp == 'IPC_protein': pI_tab.append(pI)
                
                # Patrickios model, is highly simplified so we really do not want to include this in average
                # Note: avg pI in publication (e.g., Table 1) is calculated without IPC models, for the comparison of IPC to the rest of pKa sets
                #       here avg pI is already with IPC models, yet average is ... average so should be not used so much
                if input_pKa_set_tmp != 'Patrickios': avg_ip_tab.append(pI)
                    
                ip_tab_tmp.append(pI)
            avg_ip = sum(avg_ip_tab)/len(avg_ip_tab)
            ip_string = ', '.join([str(round(n,3)) for n in ip_tab_tmp])+', '+str(round(avg_ip, 3))
            
            new_fasta_with_pI += header+os.linesep+ip_string+os.linesep+sequence+os.linesep
    
    #try to make heat map for plot
    #print(plot_libraries)
    if plot_libraries:
        try: make_heat_map(mw_tab, pI_tab, plot_file, input_pKa_set)
        except: pass
        
    print('Output file: '+output_file+os.linesep)
    f = open(output_file, 'w')
    f.write(new_fasta_with_pI)
    f.close()
   
