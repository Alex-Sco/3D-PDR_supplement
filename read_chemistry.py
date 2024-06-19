#read_chemistry.py
# This script reads the formation and destruction reactions in *.chemistry.fin files. It links them with certain colors, and plots the contribution rates of reactions (i.e. 'rate' for each reaction in the *.chemistry.fin files) along with Av values. 
# Yichen Sun yichen_sun@smail.nju.edu.cn 2024.6.17


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys

#write a example
#make a small video to show how it run


## Initial settings  
COLORls=["red","blue","purple","magenta","darkgreen","deepskyblue","limegreen","olive","black","rosybrown","navy","darkcyan","darkgoldenrod","palevioletred","chocolate","slategray","steelblue","turquoise","darkorange","slateblue","lightseagreen"]
#the colors used for plot reactions
#currenly limited to 21, should add more if there is more than 21 reactions for one species.

colorpwd="./chemical_colors/"   #the directory where the color-reactions file is stored.
modelpwd="./" #the directory where the model results is stored. 


#The names of the three models to be compared
species_for_plot=sys.argv[1]
modelname1=sys.argv[2]
modelname2=sys.argv[3]
modelname3=sys.argv[4]
figname=sys.argv[5]



# Fuction merge_reaction:
# In the first part of the Function read_combine, the texts of an equation are 'broken up', this function merge them back. 

def merge_reaction(reaction_line):
    temp=''
    for i in reaction_line:
        if i[:4] != "Rate":
            temp=temp+" "+i
        else:
            break
    return temp[1:]



# Function read_combine
# This function extracts and arranges each parameter, including the reactions and their rates.

def read_combine(mol_combine):

    Gridpoint=[]
    Av=[]
    t=[]
    nH=[]
    Tgas=[]
    abundance=[]
    Forrate=[]  # formation rate
    Desrate=[]  # destruction rate
    For_reactions=[]
    Des_reactions=[]


    for i in range(len(mol_combine)):  
        Gridpoint.append(mol_combine[i][0][2])
        Av.append(eval(mol_combine[i][0][-2]))
        t.append(eval(mol_combine[i][0][5]))
        nH.append(eval(mol_combine[i][1][2]))
        Tgas.append(eval(mol_combine[i][1][-2]))
        abundance.append(eval(mol_combine[i][-3][-1]))
        Forrate.append(eval(mol_combine[i][-2][-2]))
        Desrate.append(eval(mol_combine[i][-1][-2]))

        For_reactions_once=[]
        Des_reactions_once=[]


        count_for=0
        count_des=0

        for j in range(len(mol_combine[i])):

            if "-->" in mol_combine[i][j]:
                if mol_combine[i][j][-1]=="Rate:-100%":
                    count_des=count_des+1
                    Des_reactions_once.append([merge_reaction(mol_combine[i][j]),-100])
                elif float(mol_combine[i][j][-1][:-1])<0:
                    count_des=count_des+1
                    Des_reactions_once.append([merge_reaction(mol_combine[i][j]),float(mol_combine[i][j][-1][:-1])])
                elif float(mol_combine[i][j][-1][:-1])>0:
                    count_for=count_for+1
                    For_reactions_once.append([merge_reaction(mol_combine[i][j]),float(mol_combine[i][j][-1][:-1])])


        #sometimes, there will be no reaction output as the formation/destruction reactions in the chemistry.fin files. if that happens, label as 'no information'

        if count_for==0:
            For_reactions.append([["No Information",100]])
        else:
            For_reactions.append(For_reactions_once)

        if count_des==0:
            Des_reactions.append([["No Information",-100]])
        else:
            Des_reactions.append(Des_reactions_once)

        

    return Gridpoint,Av,t,nH,Tgas,abundance,Forrate,Desrate,For_reactions,Des_reactions

def transpose_list(matrix):  # use zip to convert row and columns, ref: https://blog.51cto.com/u_16175522/7466174
    return [list(row) for row in zip(*matrix)]



# Function arrange_reactions
#   This function arranges and combines the information of reactions and rates among all Av ranges for the plot
#   This function checks all the reactions that happened in the model among all Av ranges, if the reaction does not happen at a certain Av value, set its rate as zero.

def arrange_reactions(Av,For_reaction,Des_reaction):
    For_reaction_ls=[]
    Des_reaction_ls=[]
    For_reaction_rate_ls=[]
    Des_reaction_rate_ls=[]
    for i in range(len(Av)):
        for reac in For_reaction[i]:
            if reac[0] in For_reaction_ls:
                continue
            else: 
                For_reaction_ls.append(reac[0])
        for reac in Des_reaction[i]:
            if reac[0] in Des_reaction_ls:
                continue
            else:
                Des_reaction_ls.append(reac[0])
    
    print(For_reaction_ls)
    print(Des_reaction_ls)

    for k in range(len(For_reaction_ls)):
        For_reaction_rate_ls.append([])


    for k in range(len(Des_reaction_ls)):
        Des_reaction_rate_ls.append([])

    

    for i in range(len(Av)):
        for k in range(len(For_reaction_ls)):
            temp=transpose_list(For_reaction[i])
            if For_reaction_ls[k] in temp[0]:
                item=temp[0].index(For_reaction_ls[k])
                For_reaction_rate_ls[k].append(temp[1][item])
            else:
                For_reaction_rate_ls[k].append(0)
                

        for k in range(len(Des_reaction_ls)):
            temp=transpose_list(Des_reaction[i])
            if Des_reaction_ls[k] in temp[0]:
                item=temp[0].index(Des_reaction_ls[k])
                Des_reaction_rate_ls[k].append(temp[1][item])
            else:
                Des_reaction_rate_ls[k].append(0)

    return For_reaction_ls,For_reaction_rate_ls,Des_reaction_ls,Des_reaction_rate_ls

# Function chemistry:
#   arrange the contents from the chemistry.fin files and output the useful parameters.
#   The Gridplot, Av, t, nH, Tgas, abundance, Forrate, and Desrate are one-dimension lists, corresponding to the function of these values along the Av range
#   the For_reaction_ls and Des_reaction_ls, as one-dimension lists, storing the formation/destruction reactions happen for this species among all Av range
#   the For_reaction_rate_ls and Des_reaction_rate_ls are two-dimension lists. They have the number of items the same as the number of reactions, each item is a function of reaction rates along the Av range


def chemistry(filename,species_name):

    fp=open(filename,"r")
    ls=fp.readlines()   # read the content in the file by lines
    fp.close()
    

    #exclude the blank lines
    lssp=[]
    for i in range(len(ls)):
        a=ls[i].strip("\n").split()
        if a != []:
            lssp.append(a)
    
    
   
    species_combine=[]
    #'species_combine' instore all the useful informations for one species, every Av value, in *.chemistry.fin file


    for i in range(len(lssp)):
        if lssp[i][0]=="Species":
            if lssp[i][2]==species_name:
                count=0
                species_temp=[]
                species_temp.append(lssp[i-2])
                species_temp.append(lssp[i-1])
                for t in range(50):   # here the number 50 is an estimate of the upper limit of total lines to show reactions and related values for one Av value, one species, in the chemistry.fin file. If the reactions is too much, maybe it is not enough
                    if lssp[i+t] !=['--------------------------------------------------------------------------------'] and lssp[i+t] != []:
                        species_temp.append(lssp[i+t])
                    elif lssp[i+t] == ['--------------------------------------------------------------------------------']:
                        break
                species_combine.append(species_temp)


    #print(species_combine[0])  #if free this command we can check the structure of species_combine

   
    Gridpoint,Av,t,nH,Tgas,abundance,Forrate,Desrate,For_reactions,Des_reactions=read_combine(species_combine)
    For_reaction_ls,For_reaction_rate_ls,Des_reaction_ls,Des_reaction_rate_ls=arrange_reactions(Av,For_reactions,Des_reactions)


    return Gridpoint,Av,t,nH,Tgas,abundance,Forrate,Desrate,For_reaction_ls,For_reaction_rate_ls,Des_reaction_ls,Des_reaction_rate_ls


# Function reaction_merge
#   For different model results of the same species, merge the reactions we get from files together

def reaction_merge(ls1,ls2,ls3):
    merge=[]
    for i in ls1:
        merge.append(i)

    for i in ls2:
        if i in merge:
            continue
        else:
            merge.append(i)
    for i in ls3:
        if i in merge:
            continue
        else:
            merge.append(i)
    return merge



# Function readplotchemistry
#   This function reads the chemistry reactions, links them with colors for plot (for each species), and plots them. 
#   species_name: the name of the species for which you want to plot the formation and destruction reactions (e.g., 'CO'), should be consistent with those in the '.chemistry.fin' files.
#   model1, model2, model3: the names of three models you want to compare (the name is defined as the text before '.chemistry.fin')
#   figname: the labels of the output figures. The total figname will be [species_name]_[figname].pdf

def readplotchemistry(species_name,model1,model2,model3,figname):
    Gridpoint,Av,t,nH,Tgas,abundance,Forrate,Desrate,For_reaction_ls1,For_reaction_rate_ls1,Des_reaction_ls1,Des_reaction_rate_ls1=chemistry(modelpwd+model1+".chemistry.fin",species_name)
    Gridpoint_C1,Av_C1,t_C1,nH_C1,Tgas_C1,abundance1_C1,Forrate1_C1,Desrate1_C1,For_reaction_ls1_C1,For_reaction_rate_ls1_C1,Des_reaction_ls1_C1,Des_reaction_rate_ls1_C1=chemistry(modelpwd+model2+".chemistry.fin",species_name)
    Gridpoint_Cup,Av_Cup,t_Cup,nH_Cup,Tgas_Cup,abundance1_Cup,Forrate1_Cup,Desrate1_Cup,For_reaction_ls1_Cup,For_reaction_rate_ls1_Cup,Des_reaction_ls1_Cup,Des_reaction_rate_ls1_Cup=chemistry(modelpwd+model3+".chemistry.fin",species_name)

    print("!!!!!----",For_reaction_rate_ls1[0:2])

    For_reaction_ls1_merge=reaction_merge(For_reaction_ls1,For_reaction_ls1_C1,For_reaction_ls1_Cup)
    Des_reaction_ls1_merge=reaction_merge(Des_reaction_ls1,Des_reaction_ls1_C1,Des_reaction_ls1_Cup)
    
    # prepare a one-to-one reaction-color list and store them into files. Each color-reaction file for one species
    # output reactions to make sure all the reactions from another .py file share the same color for one species 
    # if the file is not exist, create one.

    print("Species:"+species_name,"Colorfile exists?",os.path.exists(colorpwd+"For_reaction_colors_"+species_name+".dat"))
    
    if os.path.exists(colorpwd+"For_reaction_colors_"+species_name+".dat")==False:
    
        doc=open(colorpwd+"For_reaction_colors_"+species_name+".dat","w")
        print("reactions,colors",file=doc)
        for i in range(len(For_reaction_ls1_merge)):
            print("{},{}".format(For_reaction_ls1_merge[i],COLORls[i]),file=doc)
        doc.close()
        
        doc=open(colorpwd+"Des_reaction_colors_"+species_name+".dat","w")
        print("reactions,colors",file=doc)
        for i in range(len(Des_reaction_ls1_merge)):
            print("{},{}".format(Des_reaction_ls1_merge[i],COLORls[i]),file=doc)
        doc.close()


    # --add reactions in the same file--
    
    data=pd.read_csv(colorpwd+"For_reaction_colors_"+species_name+".dat",sep=",")
    
    For_reaction_readls1=list(data["reactions"])
    colorls1=list(data["colors"])

    count=0
    for i in range(len(For_reaction_ls1_merge)):
        if For_reaction_ls1_merge[i] in For_reaction_readls1:
            continue
        else:
            count=count+1
            #print("count=",count)
            doc=open(colorpwd+"For_reaction_colors_"+species_name+".dat","a")
            print("{},{}".format(For_reaction_ls1_merge[i],COLORls[len(For_reaction_readls1)+count-1]),file=doc)
            doc.close()
    
    data=pd.read_csv(colorpwd+"Des_reaction_colors_"+species_name+".dat",sep=",")
    
    Des_reaction_readls1=list(data["reactions"])
    colorls1=list(data["colors"])

    count=0
    for i in range(len(Des_reaction_ls1_merge)):
        if Des_reaction_ls1_merge[i] in Des_reaction_readls1:
            continue
        else:
            count=count+1
            doc=open(colorpwd+"Des_reaction_colors_"+species_name+".dat","a")
            print("{},{}".format(Des_reaction_ls1_merge[i],COLORls[len(Des_reaction_readls1)+count-1]),file=doc)
            doc.close()
        
    
    # --read again, and apply color in the plots------
    
    data=pd.read_csv(colorpwd+"For_reaction_colors_"+species_name+".dat",sep=",")
    
    For_reaction_all_ls1=list(data["reactions"])
    For_reaction_all_color_ls1=list(data["colors"])
    
    
    data=pd.read_csv(colorpwd+"Des_reaction_colors_"+species_name+".dat",sep=",")
    
    Des_reaction_all_ls1=list(data["reactions"])
    Des_reaction_all_color_ls1=list(data["colors"])
    
    
    
    # ##### Below this line is the part for plot, which can be modified a lot to fit the requirements ##### 
    # --plot reactions---
    
    fig,axes=plt.subplots(2,3,figsize=(28,18),sharex=True,sharey=True)
    
    plt.subplots_adjust(wspace=0,hspace=0)
    ax1=axes[0][0]
    
    
    for i in range(len(For_reaction_ls1)):
        colorset=For_reaction_all_color_ls1[For_reaction_all_ls1.index(For_reaction_ls1[i])]
        ax1.plot(Av,For_reaction_rate_ls1[i],label=For_reaction_ls1[i],color=colorset,lw=2.5)
    
    ax1.legend(fontsize=14)
    
    ax1.set_xlabel("Av (mag)",fontsize=18)
    ax1.set_ylabel("Formation Reaction Percentage",fontsize=18)
    
    ax1.tick_params(which='minor',length=3,width=2,color="k")
    ax1.tick_params(which='major',length=8,width=2,color="k")
    
    ax1.tick_params(labelsize=18)
    ax1.text(0.5,1.05,model1,fontsize=18,ha="center",va="center",transform=ax1.transAxes)
    
    ax2=axes[0][1]
    
    for i in range(len(For_reaction_ls1_C1)):
        colorset=For_reaction_all_color_ls1[For_reaction_all_ls1.index(For_reaction_ls1_C1[i])]
        ax2.plot(Av,For_reaction_rate_ls1_C1[i],label=For_reaction_ls1_C1[i],color=colorset,lw=2.5)
    
    ax2.legend(fontsize=14)
    
    ax2.set_xlabel("Av (mag)",fontsize=18)
    ax2.set_xscale("log")
    
    ax2.tick_params(which='minor',length=3,width=2,color="k")
    ax2.tick_params(which='major',length=8,width=2,color="k")
    
    ax2.tick_params(labelsize=18)
    ax2.text(0.5,1.05,model2,fontsize=18,ha="center",va="center",transform=ax2.transAxes)
    

    # Here, by free the following command, can add and modify the annotations in the figure.
    # I use this command just to label the cosmic ray ionization rate, and the density on the figure.

    #ax2.text(0.5,1.12,r"$\rm \zeta=10^{-"+CR+r"}\ s^{-1}$  Density=$10^{3}\ \rm cm^{-3}$ "+species_name,fontsize=20,ha="center",va="center",transform=ax2.transAxes)
    
    
    
    ax3=axes[0][2]
    
    for i in range(len(For_reaction_ls1_Cup)):
        colorset=For_reaction_all_color_ls1[For_reaction_all_ls1.index(For_reaction_ls1_Cup[i])]
        ax3.plot(Av,For_reaction_rate_ls1_Cup[i],label=For_reaction_ls1_Cup[i],color=colorset,lw=2.5)
    
    ax3.legend(fontsize=14)
    
    ax3.set_xlabel("Av (mag)",fontsize=18)
    ax3.set_xscale("log")
    
    ax3.tick_params(which='minor',length=3,width=2,color="k")
    ax3.tick_params(which='major',length=8,width=2,color="k")
    
    ax3.tick_params(labelsize=18)
    
    ax3.text(0.5,1.05,model3,fontsize=18,ha="center",va="center",transform=ax3.transAxes)
    
    
    
    ax4=axes[1][0]
    
    
    for i in range(len(Des_reaction_ls1)):
        print(Des_reaction_ls1[i])
        colorset=Des_reaction_all_color_ls1[Des_reaction_all_ls1.index(Des_reaction_ls1[i])]
        ax4.plot(Av,-1*np.array(Des_reaction_rate_ls1[i]),label=Des_reaction_ls1[i],color=colorset,lw=2.5)
    
    ax4.legend(fontsize=14)
    
    ax4.set_xlabel("Av (mag)",fontsize=18)
    ax4.set_ylabel("Destruction Reaction Percentage",fontsize=18)
    
    ax4.set_xscale("log")
    
    ax4.tick_params(which='minor',length=3,width=2,color="k")
    ax4.tick_params(which='major',length=8,width=2,color="k")
    
    ax4.tick_params(labelsize=18)
    
    
    ax5=axes[1][1]
    
    for i in range(len(Des_reaction_ls1_C1)):
        colorset=Des_reaction_all_color_ls1[Des_reaction_all_ls1.index(Des_reaction_ls1_C1[i])]
        ax5.plot(Av,-1*np.array(Des_reaction_rate_ls1_C1[i]),label=Des_reaction_ls1_C1[i],color=colorset,lw=2.5)
    
    ax5.legend(fontsize=14)
    
    ax5.set_xlabel("Av (mag)",fontsize=18)
    ax5.set_xscale("log")
    
    ax5.tick_params(which='minor',length=3,width=2,color="k")
    ax5.tick_params(which='major',length=8,width=2,color="k")
    
    ax5.tick_params(labelsize=18)
    
    
    ax6=axes[1][2]
    
    for i in range(len(Des_reaction_ls1_Cup)):
        colorset=Des_reaction_all_color_ls1[Des_reaction_all_ls1.index(Des_reaction_ls1_Cup[i])]
        ax6.plot(Av,-1*np.array(Des_reaction_rate_ls1_Cup[i]),label=Des_reaction_ls1_Cup[i],color=colorset,lw=2.5)
    
    ax6.legend(fontsize=14)
    
    ax6.set_xlabel("Av (mag)",fontsize=18)
    ax6.set_xscale("log")
    
    ax6.tick_params(which='minor',length=3,width=2,color="k")
    ax6.tick_params(which='major',length=8,width=2,color="k")
    
    ax6.tick_params(labelsize=18)
    
    #save the figure with a name
    plt.savefig(species_name+"_"+figname+".pdf")
    plt.close()
    


##--- The main function --- plot the figures, this part can be modified to meet the requirements.

readplotchemistry(species_for_plot,modelname1,modelname2,modelname3,figname)


#here is also an example for batchly use the readplotchemistry function:

'''
CRls=[15,17]
species_ls=["CO","OH"]#,"OH","HCO+","H3+","H2O","H2O+","CH","OH+"]

for i in CRls:
    for j in species_ls:
        CRrate=i
        CR=str(CRrate)
        readplotchemistry(j,"O0C0D0_"+CR,"O0C1D0_"+CR,"O0CupD0_"+CR,"O0models")
'''

#merge pdfs, use 'pdfunite'
#
#for j in species_ls:
#    os.system("pdfunite "+j+"_*.pdf all_"+j+"_chemical_reactions_test.pdf")





