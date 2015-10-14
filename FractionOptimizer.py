import sys
import matplotlib
matplotlib.use('Agg')
import csv
import random
import pylab
import seaborn

from datetime import timedelta
from scipy import stats
from pyteomics import biolccc, mass
from pyteomics import pylab_aux
import numpy as np

import ConfigParser
Config = ConfigParser.ConfigParser()
Config.read(sys.argv[1])
Config.sections()

#function for prediction the retation time of peptides in fractionation by BioLCCC
def biolccc_func():
    frac=[]
    res=[]
    pep=[]
    RT_res=[]
    dic_theor={}
    dic_exp={}
    peptides_theor=[]
    ColumnLength=float(Config.get('BioLCCC', 'ColumnLength'))
    ColumnDiameter=float(Config.get('BioLCCC', 'ColumnDiameter'))
    ColumnPoreSize=float(Config.get('BioLCCC', 'ColumnPoreSize'))
    SecondSolventConcentrationA=float(Config.get('BioLCCC', 'SecondSolventConcentrationA'))
    SecondSolventConcentrationB=float(Config.get('BioLCCC', 'SecondSolventConcentrationB'))
    FlowRate=float(Config.get('BioLCCC', 'FlowRate'))
    ColumnRelativeStrength=float(Config.get('BioLCCC', 'ColumnRelativeStrength'))
    time_grad_fr_str=Config.get('BioLCCC', 'time_grad_fr')
    time_grad_fr=[float(x) for x in time_grad_fr_str.split(',')]
    percent_grad_fr_str=Config.get('BioLCCC', 'percent_grad_fr')
    percent_grad_fr=[float(x) for x in percent_grad_fr_str.split(',')]
    with open(sys.argv[2], 'rb') as fIn: 
        reader = csv.reader(fIn, delimiter='\t')
        next(reader)
        for row in reader:
            frac.append(row[0])
            dic_exp[row[0]]=row[6]
    a=len(frac)
    myChromoConditions = biolccc.ChromoConditions()
    myChromoConditions.setColumnRelativeStrength(ColumnRelativeStrength)
    myChromoConditions.setColumnLength(ColumnLength)
    myChromoConditions.setColumnDiameter(ColumnDiameter)
    myChromoConditions.setColumnPoreSize(ColumnPoreSize)
    myChromoConditions.setSecondSolventConcentrationA(SecondSolventConcentrationA)
    myChromoConditions.setSecondSolventConcentrationB(SecondSolventConcentrationB)
    myChromoConditions.setFlowRate(FlowRate)


    myGradient = biolccc.Gradient()
    for i in range(len(time_grad_fr)):
        myGradient.addPoint(time_grad_fr[i], percent_grad_fr[i])
    myChromoConditions.setGradient(myGradient)


    standard_aminoacids = set(k for k in mass.std_aa_comp if '-' not in k)

    for peptide in frac:
        try:
            RT = biolccc.calculateRT(peptide,
                biolccc.rpAcnTfaChain,
                myChromoConditions,0,True,False)
            RT_res.append(RT)
            peptides_theor.append(peptide)
            dic_theor[peptides]=RT
        except:
            pass
    return (dic_exp,frac,RT_res,peptides_theor,dic_theor,myChromoConditions)


#function for calculate time of nonsorb. component by BioLCCC
def nonsorb():
    nonsorb=[]
    nonsorb.append('Ac-O-NH2')
    myChemicalBasis = biolccc.ChemicalBasis(
        biolccc.RP_ACN_FA_ROD)
    myChemicalBasis.addChemicalGroup(
        biolccc.ChemicalGroup(
       'Nonsorb',      
        'O',                 
       0.0,                  
       0,      
       0))  
    for peptide in nonsorb:
            time_nonsorb_biolccc = biolccc.calculateRT(peptide,
                myChemicalBasis,
                myChromoConditions,0,True,False)
    return time_nonsorb_biolccc


#test of the peptides distribution in optim. fractionation
def test(RT_res):
 #   print stats.scoreatpercentile(RT_res, [20,80])
    pylab.hist(RT_res, bins=np.arange(0,40,1),normed=0,alpha=0.3,facecolor='blue')
    pylab.xlim(0,40)
    pylab.xlabel('time,min')
    pylab.ylabel('amount of pep')
    pylab.title('Gradient for fractionation')
    pylab.savefig('%s/test_optimal_gradient.jpg'%sys.argv[3])


#optimizer
def optimizer():
    time_fraction=[]
    optim_fraction=[]
    amount_peptides=[]
    delta_optim=[]
    time=[]
    number_frac=[]
    dic_number_fr={}
    delay_time_grad=float(Config.get('General', 'delay time grad'))
    fr=int(Config.get('General', 'Amount of fractions'))
    time_nonsorb_exp=float(Config.get('General', 'time_nonsorb_exp'))
    time_det_collector=float(Config.get('General', 'time_det_collector'))
    pylab.figure(figsize=(6,3))
    percentile=np.linspace(0, 100, fr+1)
    time=stats.scoreatpercentile(RT_res, percentile[1:])
    time_fraction=list(time)
    for i in range(len(time_fraction)):
        time_fraction[i]=round(time_fraction[i],1)
    time_fraction.insert(0,0)
    full_time=float(Config.get('General', 'Full gradient time for fractionation'))
    time_fraction[fr]=float(full_time)
    for j in range(fr):
        delta_optim.append(round((time_fraction[j+1]-time_fraction[j]),1))
    number_frac=np.digitize(RT_res, time_fraction, right=False)
    dic_number_fr = dict(zip(peptides_theor, number_frac))
    #correction (delay time,nonsorb.component, etc.)
    time_fraction_change=[]
    pylab.rcParams['font.family'] = 'DejaVu Sans'
    events, edges, patches = pylab.hist(RT_res,time_fraction,normed=0,facecolor='blue', alpha=0.5)
    pylab.xlabel('t,min')
    pylab.ylabel('amount of peptides')
    pylab.title('Time for collection of the fractions')
    pylab.savefig('%s//optimal_gradient.jpg'%sys.argv[3])
    for i in range(len(time_fraction)):
        time_fraction_change.append(round(time_fraction[i]+delay_time_grad-time_nonsorb_exp-nonsorb()+time_det_collector,3))
    time_fraction_change[0]='time of nonsorb'
    time_fraction_change.pop()
    f=open('%s/optimal_fractions.txt'%sys.argv[3],'w')
    print >>f, time_fraction_change, 'start time for collection of the fraction'
    print >>f, delta_optim, 'time for collection of the fraction (delta)'
    print >>f, events, 'amount of pep in each fraction'
    f.close()
    return time_fraction



#peptides distribution in the fractions
def pep_distr(RT_res,time_fraction):
    fr=int(Config.get('General', 'Amount of fractions'))
    delay_time_analytical=float(Config.get('General', 'delay_time_analytical'))
    percent_acn=float(Config.get('General', 'percent_acn'))
    RT=[]
    RT_fr=[]
    number_frac=[]
    peptides=[]
    percent_ACN=[]
    peptides=np.array(frac)
    RT=np.array(RT_res)
    number_frac=np.digitize(RT, time_fraction, right=False)
    pylab.figure(figsize=(20,10))
    for i in range(fr):
        RT_fr=[]
        for pep in peptides[number_frac==i+1]:
            RT_fr.append(float(dic_exp[str(pep)]))
        pylab.subplot(2,3,i+1)    
        pylab.hist(RT_fr, bins=np.arange(0,300,10),normed=0,alpha=0.3,facecolor='blue')
        pylab.title('fraction%d'%(i+1))
        pylab.xlim(0,300)
        t1 = percent_acn*grad(stats.scoreatpercentile(RT_fr, [5])[0]-delay_time_analytical)
        t2 = percent_acn*grad(stats.scoreatpercentile(RT_fr, [95])[0]-delay_time_analytical)
        percent_ACN.append((t1,t2))
        pylab.savefig('%s/pep_distribution.jpg'%sys.argv[3]) 
    return percent_ACN


#additional functions
def lin_func(x1,y1,x2,y2,x): #return equation of each linear parts of the gradient
    a=(y1-y2)/(x1-x2)
    b=y1-a*x1
    y=a*float(x)+b
    return float(y)
def grad(time): #return percent of ACN for x min 
    time_grad_str=Config.get('General', 'time_grad')
    time_grad=[float(x) for x in time_grad_str.split(',')]
    percent_grad_str=Config.get('General', 'percent_grad')
    percent_grad=[float(x) for x in percent_grad_str.split(',')]
    if time>0:
        num=np.digitize([time], time_grad, right=False)
        num_piece=num[0]
        x1=time_grad[num_piece-1]
        x2=time_grad[num_piece]
        y1=percent_grad[num_piece-1]
        y2=percent_grad[num_piece]
        return lin_func(x1,y1,x2,y2,time)
    else:
        return float(percent_grad[0])


#analitical gradient for fractions
def grad_for_frac(error):
    time_str=Config.get('General', 'time_str')
    time=[float(x) for x in time_str.split(',')]
    fr=int(Config.get('General', 'Amount of fractions'))
    pylab.figure(figsize=(20,13))
   
    percent_ACN=pep_distr(RT_res,time_fraction)
    f=open('%s/analytical_grad_for_fractions.txt'%sys.argv[3],'w')
    pylab.figure(figsize=(20,10))
    for i in range(fr):
        acn_f=[]
        acn=[3,3,70,70,3,3]
        if percent_ACN[i][0]>=acn[0]:
            acn.insert(2,round(percent_ACN[i][0],3))
        else:
            acn.insert(2,acn[0])
        acn.insert(3,round(percent_ACN[i][1]+error,3))
        acn_f.append(acn)
        print >>f, acn, 'fraction%d'%(i+1)
        pylab.subplot(2,3,i+1)
        pylab.plot(time,acn,c='blue',linewidth=3.0,alpha=0.5)
        pylab.ylim(0,100)
        pylab.title('fraction%d'%(i+1))
        pylab.xlabel("time,min")
        pylab.ylabel('percent_ACN')
        pylab.savefig('%s/analytical_grad_for_fractions.jpg'%sys.argv[3])


if __name__ == '__main__':
    dic_exp,frac,RT_res,peptides_theor,dic_theor,myChromoConditions=biolccc_func()
    time_nonsorb_biolccc=nonsorb()
    test(RT_res)
    time_fraction=optimizer()
    percent_ACN=pep_distr(RT_res,time_fraction)
    grad_for_frac(10)
    print ('Enjoy your optimal fractions!')
