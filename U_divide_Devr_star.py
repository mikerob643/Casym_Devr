# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 10:54:31 2023

@author: mjrobinson
"""


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


import statsmodels.api as sm

from scipy.stats import binned_statistic


plt.rcParams.update({'font.size': 9})
plt.rcParams['font.family'] = ['serif']
plt.rcParams['font.serif'] = ['Times New Roman']
# this is not needed just for plotting preference
from cmcrameri import cm



df_casym= pd.read_csv('Casym2Mikey_ForEachRidge4.csv',header=0)
#data from matlab code 


caysm_dist=np.array((df_casym.distance))


caysm_median=np.array((df_casym.Casym_medians))
caysm_error=np.array((df_casym.Casym_std/np.sqrt(df_casym.numProfiles)))
# -1.09455478666666672 is the average value for the last 6 ridges in dragons back 


caysm_median=caysm_median-0.09455478666666672

caysm_number=np.array((df_casym.numProfiles))

# digitized hillslope length from Hurst DBPR paper 
df_hurst_length= pd.read_csv('hurst_hillslope_length.csv',header=0)
xl=np.array((df_hurst_length.x))

yl=np.array((df_hurst_length.hillslope_length))









df= pd.read_csv('DragonsBack_Ridgelines_1m_Corrected.csv',header=0,float_precision='round_trip')


basin_data = df

basin_data.columns = ['x', 'y', 'z','ridge']

basin_data = basin_data.dropna()
xridge = np.array(df['x'])
yridge = np.array(df['y'])
zridge = np.array(df['z'])
ridge  =np.array(df["ridge"])
x1ridge = np.array(df['x'])
y1ridge = np.array(df['y'])
z1ridge = np.array(df['z'])


#%%
figure_proportions = (9, 6)
plt.figure(figsize=figure_proportions)
uniqueridges=np.unique(ridge)
# plt.scatter(xridge,yridge)
for i in range (np.size(np.unique(ridge))):
    print(i)
    plt.scatter(xridge[ridge==uniqueridges[i]],yridge[ridge==uniqueridges[i]],label="{}".format(uniqueridges[i]))
    plt.legend()
#%%
basin_list = basin_data.ridge.unique()


x1=xridge[ridge==np.max(ridge)]
y1=yridge[ridge==np.max(ridge)]


x1=x1[y1==np.min(y1)]
y1=y1[y1==np.min(y1)]


plt.scatter(xridge,yridge)
plt.scatter(x1,y1)
#%%

whatwecut2=np.array((6,8,10,11,12,13,15,17,18,19,20,22,23,25,27,29,31,32,33,42,44,45,46,48,50,52,54,56,58,74,75))
# these are divides that do not meet our criteria 
whatwecut2=np.sort(whatwecut2)
for i in range(np.size(whatwecut2)):
    basin_list=basin_list[basin_list!=whatwecut2[i]]
#%%

#find distance 

plt.figure(figsize=(10,10))  # Adjust figure size as needed

for i in range (np.size(basin_list)):
    print(i)
    plt.scatter(xridge[ridge==basin_list[i]],yridge[ridge==basin_list[i]],label="{}".format(basin_list[i]))
    plt.legend()
plt.savefig('check.jpg',dpi=1000,bbox_inches='tight')
# look at divides 
#%%

for i in range(np.size(whatwecut2)):
    basin_list=basin_list[basin_list!=whatwecut2[i]]

#%%


area=np.zeros(np.size(basin_list))

normal_area=np.zeros(np.size(basin_list))

left=np.zeros(np.size(basin_list))
right=np.zeros(np.size(basin_list))

normleft=np.zeros(np.size(basin_list))
normright=np.zeros(np.size(basin_list))

leftavg=np.zeros(np.size(basin_list))
rightavg=np.zeros(np.size(basin_list))


xisforgettingl=np.zeros(np.size(basin_list))

distance=np.zeros(np.size(basin_list))


ridgexpos=np.zeros(np.size(basin_list))
ridgeypos=np.zeros(np.size(basin_list))
ridgexneg=np.zeros(np.size(basin_list))
ridgeyneg=np.zeros(np.size(basin_list))

kristin_x=()
kristin_y=()
kristin_dist=()
kristin_ridge_id=()

for f in range(0, np.size(basin_list), 1):


    basinpick = basin_list[f]
    grp2=df[(df["ridge"]==basin_list[f])]
    

    xridge = np.array(grp2['x'])
    yridge = np.array(grp2['y'])
    
    if basinpick==27:
        xridge = xridge[:-10]
        yridge = yridge[:-10]

    


    distance[f]=+np.sqrt(((x1[0])-(xridge[yridge==np.min(yridge)][0]))**2
                         +((y1[0])-(yridge[yridge==np.min(yridge)][0]))**2)



    kristin_x=np.append(kristin_x,xridge)
    kristin_y=np.append(kristin_y,yridge)
    kristin_dist=np.append(kristin_dist,np.full(np.size(xridge),distance[f]))
    kristin_ridge_id=np.append(kristin_ridge_id,np.full(np.size(xridge),basinpick))




#%%

# # export all of these values so I dont have to run this above cell again. it takes awhile 
k_data=pd.DataFrame()
k_data["time"]= np.zeros(np.size(kristin_dist))

k_data["distance"]= kristin_dist
k_data["ridge"]= kristin_ridge_id

k_data["x"]=      kristin_x
k_data["y"]= kristin_y

k_data["z"]=  np.zeros(np.size(kristin_dist))


k_data.to_csv("k_data.csv")
  

# this can be ignored and was only used for our specific case 

#%%
for f in range(0, np.size(basin_list), 1):


    basinpick = basin_list[f]
    grp2=df[(df["ridge"]==basin_list[f])]
    

    xridge = np.array(grp2['x'])
    yridge = np.array(grp2['y'])


    distance[f]=+np.sqrt(((x1[0])-(xridge[yridge==np.min(yridge)][0]))**2
                          +((y1[0])-(yridge[yridge==np.min(yridge)][0]))**2)


    findtopandbottom=xridge+yridge
    figure_proportions = (9, 6)
    plt.figure(figsize=figure_proportions)
    fi, ax = plt.subplots()
    plt.title(basinpick)
    
    bottomx=xridge[findtopandbottom==np.min(findtopandbottom)][0]
    bottomy=yridge[findtopandbottom==np.min(findtopandbottom)][0]
    
    
    topx=xridge[findtopandbottom==np.max(findtopandbottom)][0]
    topy=yridge[findtopandbottom==np.max(findtopandbottom)][0]
     
    
    plt.scatter(topx,topy,color="black",zorder=11,s=40)
    plt.scatter(bottomx,bottomy,color="red",zorder=10,s=40)
    plt.plot(xridge,yridge)
    
    normaldist=np.sqrt((topx-bottomx)**2+(topy-bottomy)**2)
    
 
    
 
    
    m=((topy-bottomy)/(topx-bottomx))
    
    
    b=(bottomx*topy-topx*bottomy)/(bottomx-topx)
    
    newsx=np.linspace(np.min(xridge),np.max(xridge),10000 )

    
    newys=m*xridge+b
    
    
    newyss=m*newsx+b
    




    plt.scatter(xridge,newys)
    

    plt.scatter(xridge,yridge)
    
    



    add =np.zeros(np.size(newys))
    
    
    positive=np.zeros(np.size(newys))
    
    
    neg=np.zeros(np.size(newys))
    

    for g in range(np.size(newys)):
        if newys[g]>yridge[g]:
           
            findshortest=np.zeros(np.size(newyss))
            
            
            indexe=np.arange(0,np.size(newyss),1)
            
            
            for p in range(np.size(newyss)):
                
                findshortest[p]=np.sqrt((xridge[g]-newsx[p])**2+(yridge[g]-newyss[p])**2)
                
                
            indi=indexe[findshortest==np.min(findshortest)]
            indi=indi[0]
                                                    

                  
    
            positive[g]=np.min(findshortest)
            
            # in the paper this was decided to be represented as negative
        
            plt.plot((xridge[g],newsx[indi]), (yridge[g],newyss[indi]),color="red" )
            
         
            
 
            
        if newys[g]<yridge[g]:

            findshortest=np.zeros(np.size(newyss))
            
            
            indexe=np.arange(0,np.size(newyss),1)
            
            
            for p in range(np.size(newyss)):
                
                findshortest[p]=np.sqrt((xridge[g]-newsx[p])**2+(yridge[g]-newyss[p])**2)
                
                
            indi=indexe[findshortest==np.min(findshortest)]
            indi=indi[0]
                                                    

            
            neg[g]=-1*np.min(findshortest)
            
            
    
            
            plt.plot((xridge[g],newsx[indi]), (yridge[g],newyss[indi]),color="blue" )
            
  
   
        if np.size(positive[positive!=0])!=0:
                  
    
            leftavg[f]=np.median(positive)
            left[f]=np.max(positive)
            normleft[f]=np.max(positive)/normaldist
            
            
            
     

            
        if np.size(neg[neg!=0])!=0: 
         
            rightavg[f]=np.median(neg)
            
            right[f]=np.min(neg)
            normright[f]=np.min(neg)/normaldist
            
            # right and left are flipped. normright really means deviation to the left. and normleft really means deviation to the right. 
            #
    
       



        
    plt.text(.01, .99, '{}'.format(right[f]+left[f]), ha='left', va='top', transform=ax.transAxes)




# # export all of these values 
sout=pd.DataFrame()
sout["rightavg"]= rightavg

sout["right"]= right
sout["normright"]=      normright
sout["rightavg"]= rightavg

sout["left"]= left
sout["normleft"]=     normleft

sout["leftavg"]= leftavg


sout.to_csv("sout.csv")
  
    # enter in data from above loop without doint it 
   
   #%%
   
sout= pd.read_csv('sout.csv',header=0)
right= np.array((sout["right"]))
normright= np.array((sout["normright"]))
rightavg= np.array((sout["rightavg"]))

left= np.array((sout["left"]))
normleft= np.array((sout["normleft"]))

leftavg=np.array((sout["leftavg"]))


sumation=normright+normleft


casym_apples=np.zeros(np.size(leftavg))
casym_std=np.zeros(np.size(leftavg))

casym_only=np.zeros(np.size(leftavg))
casym_only_no_divide=np.zeros(np.size(leftavg))

             






        
            

#%%


# 45 gets it to hurst starting point. 

distance1=distance+45
   
time=distance1/.033




#%%



s2=pd.DataFrame()
s2["sumation"]=sumation
s2["sumation_2"]=sumation
s2["time"]=time
s2["dist"]=distance1





#%%
veldir=np.zeros(42)
veldir[:]=1
b22=np.copy(time)





save_casymss=np.zeros(np.size(caysm_median))
for i in range(np.size(time)):
    if b22[i]>60000:
        b22[i]=b22[i]-60000
        # this represents the point of transtion from growth to decay zone 

D=0.0086
lambdaa=30
tau_d = 1

# lambda using hurst hillslope length
lambdaa_new=np.zeros(np.size(distance1))
for i in range (np.size(lambdaa_new)):
    lambdaa_new[i]=np.interp(distance1[i],xl,yl)
    

t_save=np.zeros(np.size(leftavg))
u_error=np.zeros(np.size(caysm_error))
ue2=np.zeros(np.size(caysm_error))
for i in range(np.size(leftavg)):
    tdim=b22[i]
    t_save[i]=time[i]
    TD=lambdaa_new[i]**2/D
    t_hat = tdim/TD


    save_casyms=caysm_median[i]
    
    
   
    veldir=np.zeros(np.size(save_casyms))



    if  save_casyms<=1:
        veldir=-1
    if  save_casyms>1:
        veldir=1
        save_casyms=1/ save_casyms
        
        
    save_casymss[i]=save_casyms
        
             

    theta_I = save_casyms
    
    u_d = 3*(1-theta_I) / np.sqrt(6*tau_d * (2*tau_d - 3*t_hat + 2*theta_I*(3*t_hat+tau_d) + theta_I**2*(2*tau_d - 3*t_hat)))#   % Equation 36
    u_dim = u_d*(D/lambdaa_new[i])   #;  % Mudd/Furbish said lambdaa/D, but I think that is years per unit length 

    u_dim =u_dim*veldir
         
    casym_apples[i]=u_dim
 
    # casym_apples is velocitys for divide migration 

    ue2[i]=(3*tau_d*(theta_I+1))/(((theta_I**2)*(t_hat-.666667*tau_d)-2*theta_I*(t_hat+.333333*tau_d)+t_hat-.666667*tau_d)  *(np.sqrt(-18*tau_d*((theta_I**2)*(t_hat-.666667*tau_d)-2*theta_I*(t_hat+.333333*tau_d)+t_hat-.666667*tau_d))))                                                                                        
    ue2[i] = np.abs(ue2[i])*caysm_error[i]
    ue2[i]=ue2[i]*(D/lambdaa_new[i])  
    # error propegation 
    
u_dimensional= casym_apples[~np.isnan(casym_apples)]
u_error=ue2[~np.isnan(casym_apples)]

#%%


DEV_STAR=pd.DataFrame()
DEV_STAR["dev"]=sumation[~np.isnan( casym_apples)]
DEV_STAR["time"]=b22[~np.isnan( casym_apples)]
DEV_STAR["time_c"]=time[~np.isnan( casym_apples)]



dev_each=np.array((DEV_STAR.dev))

DEV_STAR_dist=np.array((distance1))


ridgeid=basin_list







