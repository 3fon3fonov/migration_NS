#!/usr/bin/python
import numpy as np
import os, sys
import matplotlib
#mpl.use('Agg')

matplotlib.use('Agg')


import matplotlib.pyplot as plt
from pylab import *

from subprocess import PIPE, Popen
#import scipy.optimize as op
#import emcee
#import signal
import corner_ES as corner
#import dynesty
import dill
dill.settings['fmode']

#import gc
import dynesty_2_0 as dynesty
#import dynesty
import dynesty_patch
dynesty.results =  dynesty_patch 


from threading import Thread
from pathos.multiprocessing import ProcessingPool as Pool

'''
Trifon Trifonov 2019
'''

THIRD = 1.0/3.0
PI    = 3.14159265358979e0
TWOPI = 2.0*PI
GMSUN = 1.32712497e20
AU=1.49597892e11
incl = 90.0
sini = np.sin(PI*(incl/180.0))
 
G  = 6.67384e-11 


########################### For nice plotting ##################################

mpl.rcParams['axes.linewidth'] = 2.0 #set the value globally
mpl.rcParams['xtick.major.pad']='1'
mpl.rcParams['ytick.major.pad']='2'


# set tick width
mpl.rcParams['xtick.major.size'] = 8
mpl.rcParams['xtick.major.width'] = 2
mpl.rcParams['xtick.minor.size'] = 5
mpl.rcParams['xtick.minor.width'] = 2

mpl.rcParams['ytick.major.size'] = 8
mpl.rcParams['ytick.major.width'] = 2
mpl.rcParams['ytick.minor.size'] = 5
mpl.rcParams['ytick.minor.width'] = 2


rc('text',usetex=True)
font = {'family' : 'normal','weight' : 'black','size': 22,'serif':['Helvetica']}
rc('font', **font)

########################## run the wrapper #####################################



class FunctionWrapper(object):
    """
    This is a hack to make the likelihood function pickleable when ``args``
    or ``kwargs`` are also included.
    """
    def __init__(self, f, args, kwargs=None):
        self.f = f
        self.args = [] if args is None else args
        self.kwargs = {} if kwargs is None else kwargs

    def __call__(self, x):
        try:
            result = self.f(x, *self.args, **self.kwargs)
            #print(x, result)
            return result
        except:  # pragma: no cover
            import traceback
            print("emcee: Exception while calling your likelihood function:")
            print("  params:", x)
            print("  args:", self.args)
            print("  kwargs:", self.kwargs)
            print("  exception:")
            traceback.print_exc()
            raise


def run_command_with_timeout(args, secs):
    proc = Popen(args, shell=True, preexec_fn=os.setsid, stdout=PIPE)
    proc_thread = Thread(target=proc.wait)
    proc_thread.start()
    proc_thread.join(secs)
    text = proc.communicate()[0]
    flag = 1
    if proc_thread.is_alive():
        print(proc.pid)
        try:
            os.killpg(proc.pid, signal.SIGTERM)
        except OSError:

            print('Process #{} killed after {} seconds'.format(proc.pid, secs))
            flag = -1
            #text = '0 0 0 0'
            return text.decode('utf-8'),flag
    #return proc, flag , text.decode('utf-8')
    return text.decode('utf-8'),flag
    
    


################  concentrate parameter in arrays ################# 
def concentrate_par(par,flag,bounds,el_str):
     
    f = [idx for idx in range(len(flag)) if flag[idx] ==1 ]

    p = []  #'p" are the fitted parameters
    b = []
    e = []

    for j in range(len(par)):
        if flag[j] > 0:
            p.append(par[j])
            b.append(bounds[j])
            e.append(el_str[j])

    return p,f,b,e 


 
 
 
################################ wrapper #######################################   

def migration(p, par,flag_ind, npl):


    if p[0] != 0:
        rand_dir = int(p[0]*np.random.randint(0,100000000000))
    else:
        rand_dir = int(p[3]*np.random.randint(0,100000000000))            
    
    os.system("mkdir %s"%str(rand_dir))        

    os.chdir("%s"%str(rand_dir))  

    for j in range(len(p)):
        par[flag_ind[j]] = p[j] 


    Stellar_mass =  1 #par[-1]
    pl_file = 'pl.in'
    mtiny =  1e-6
    t_max = 10**(abs(par[7*npl+3])) * 0.7
    dt = 0.05


    ############ Generate geninit_j.in, pl.in, param.in in Jacobi coor. #####
    ############ and start the N-body code                              #####

    getin_file = open('geninit_j.in', 'w') 

    getin_file.write("""1 
%s
%s
1.d0
%s
    """%(str(Stellar_mass),str(npl),str(pl_file)))

 
    for j in range(npl): 
        getin_file.write('%s \n'%str(par[7*j])) 
        getin_file.write('%s %s %s %s %s %s \n'%(str(par[7*j +1]),str(par[7*j +2]),str(par[7*j +5]),str(par[7*j +3]),str(par[7*j +6]),str(par[7*j +4])) ) 

    getin_file.close()

     

    ##### crate the pl.in file, (ignore the std. output). ######
    timeout_sec = 300.0
    result, flag = run_command_with_timeout("../geninit_j3 < geninit_j.in", timeout_sec)


    ##### crate the param.in file (change only the "t_max" and the "dt" for now) ######
    param_file = open('param.in', 'w') 
     
    param_file.write("""0.0d0 %s %s
3.9999999999d1 3.9999999999d4
F T T T T F
.22 20. 20. -1. T
bin.dat
unknown
    """%(str(t_max), str(dt)))
     
    param_file.close()


    ##### crate the symba_input.in file (keep it as it is for now) ######

    symba_file = open('symba_j_migrate6.in', 'w') 
    symba_file.write("""param.in
%s
%s
%s %s
-%s %s
    """%(pl_file,str(mtiny), str(par[7*npl]),str(par[7*npl+1]),str((10**(par[7*npl+3]))),str(par[7*npl+2])) )
          
    symba_file.close()

    ##### run SyMBA. ######
    timeout_sec = 86400.0 # do not run more than 1 day. up to the user!

    result, flag = run_command_with_timeout("../swift_symba5_j_migrate6 < symba_j_migrate6.in", timeout_sec)

    #os.system('../swift_symba5_j_migrate6 < symba_j_migrate6.in')


    ################## make the ascii output for the plotting ############### 


    for j in range(npl):

        os.system('''../follow_symba2 << EOF > /dev/null 2>&1
param.in
%s
%s
EOF
    '''%(pl_file, j+2))
        os.system('mv follow_symba.out follow_symba.%s.out'%(j+1))


    skip_header = 0
    skip_footer = 0

    try:
        T1    = genfromtxt('./follow_symba.1.out',skip_header=skip_header, unpack=True,skip_footer=skip_footer, usecols = [0])      
        T1 = T1/1000.0
        T2  = T1
    
    except:
        T2    = genfromtxt('./follow_symba.2.out',skip_header=skip_header, unpack=True,skip_footer=skip_footer, usecols = [0])      
        T2 = T2/1000.0
        T1  = T2    

    
    
    pl_range = 1
    range_t1 = len(T1)
    range_max = int(range_t1/pl_range) 

    omega1 =  np.genfromtxt('./follow_symba.1.out',skip_header=skip_header, unpack=True,skip_footer=skip_footer, usecols = [6])
    omega2 =  np.genfromtxt('./follow_symba.2.out',skip_header=skip_header, unpack=True,skip_footer=skip_footer, usecols = [6])

    capm1 =  np.genfromtxt('./follow_symba.1.out',skip_header=skip_header, unpack=True,skip_footer=skip_footer, usecols = [7])
    capm2 =  np.genfromtxt('./follow_symba.2.out',skip_header=skip_header, unpack=True,skip_footer=skip_footer, usecols = [7])

    a_1   = np.genfromtxt('./follow_symba.1.out',skip_header=skip_header, unpack=True,skip_footer=skip_footer, usecols = [2])
    a_2   = np.genfromtxt('./follow_symba.2.out',skip_header=skip_header, unpack=True,skip_footer=skip_footer, usecols = [2])

    e_1   = np.genfromtxt('./follow_symba.1.out',skip_header=skip_header, unpack=True,skip_footer=1, usecols = [3])
    e_2   = np.genfromtxt('./follow_symba.2.out',skip_header=skip_header, unpack=True,skip_footer=1, usecols = [3])

    per1  = omega1
    per2  = omega2

    omegai1 = np.radians(omega1)
    omegai2 = np.radians(omega2)
    capmi1 = np.radians(capm1)
    capmi2 = np.radians(capm2)


    range_max_e =range_max
    range1 = int(range_max_e)
 
    theta1max = PI
    theta1min = -PI
    theta2max = PI
    theta2min = -PI
    theta3max = PI
    theta3min = -PI
    theta4max = PI
    theta4min = -PI
    theta5max = PI
    theta5min = -PI
    theta6max = PI
    theta6min = -PI
    theta1 =np.zeros(range1)
    theta2 =np.zeros(range1)
    theta3 =np.zeros(range1)
    theta4 =np.zeros(range1)
    theta5 =np.zeros(range1)
    theta6 =np.zeros(range1)
    capomi1 = np.zeros(range1)
    capomi2 = np.zeros(range1)

    delta_theta1 =np.zeros(range1)
    delta_theta2 =np.zeros(range1)    
    delta_theta3 =np.zeros(range1)
    delta_theta4 =np.zeros(range1)
    delta_theta5 =np.zeros(range1)

    c =  np.zeros(range1)

    for f in range(0,range1):

        #theta1[f]  = (capomi1[f]  + omegai1[f]  + capmi1[f] ) + (capomi1[f]  + omegai1[f] )   - 2*(capomi2[f]  + omegai2[f]  + capmi2[f] )
        #theta2[f]  = (capomi1[f]  + omegai1[f]  + capmi1[f] ) - 2*(capomi2[f]  + omegai2[f]   + capmi2[f])     + (capomi2[f]  + omegai2[f] )

        theta1[f]  = (capomi1[f]  + omegai1[f]  + capmi1[f] )%(2.0*PI) + (5*(capomi1[f]  + omegai1[f] )%(2.0*PI))%(2.0*PI)  - (6*(capomi2[f]  + omegai2[f]  + capmi2[f] )%(2.0*PI))%(2.0*PI) 
        theta2[f]  = (capomi1[f]  + omegai1[f]  + capmi1[f] )%(2.0*PI) + (4*(capomi1[f]  + omegai1[f] )%(2.0*PI))%(2.0*PI)  - (6*(capomi2[f]  + omegai2[f]  + capmi2[f] )%(2.0*PI))%(2.0*PI) +(1*(capomi2[f] +omegai2[f] ))%(2.0*PI) 
        theta3[f]  = (capomi1[f]  + omegai1[f]  + capmi1[f] )%(2.0*PI) + (3*(capomi1[f]  + omegai1[f] )%(2.0*PI))%(2.0*PI)  - (6*(capomi2[f]  + omegai2[f]  + capmi2[f] )%(2.0*PI))%(2.0*PI) +(2*(capomi2[f] +omegai2[f] ))%(2.0*PI) 
        theta4[f]  = (capomi1[f]  + omegai1[f]  + capmi1[f] )%(2.0*PI) + (2*(capomi1[f]  + omegai1[f] )%(2.0*PI))%(2.0*PI)  - (6*(capomi2[f]  + omegai2[f]  + capmi2[f] )%(2.0*PI))%(2.0*PI) +(3*(capomi2[f] +omegai2[f] ))%(2.0*PI) 
        theta5[f]  = (capomi1[f]  + omegai1[f]  + capmi1[f] )%(2.0*PI) + (1*(capomi1[f]  + omegai1[f] )%(2.0*PI))%(2.0*PI)  - (6*(capomi2[f]  + omegai2[f]  + capmi2[f] )%(2.0*PI))%(2.0*PI) +(4*(capomi2[f] +omegai2[f] ))%(2.0*PI) 
        theta6[f]  = (capomi1[f]  + omegai1[f]  + capmi1[f] )%(2.0*PI) + (5*(capomi2[f]  + omegai2[f] )%(2.0*PI))%(2.0*PI)  - (6*(capomi2[f]  + omegai2[f]  + capmi2[f] )%(2.0*PI))%(2.0*PI) 


        theta1[f]  = degrees(theta1[f]%(2.0*PI))%360.0 
        if theta1[f] >= 180.0:
                theta1[f] = theta1[f] - 360.0
        theta2[f]  = degrees(theta2[f]%(2.0*PI))%360.0 
        if theta2[f] >= 180.0:
                theta2[f] = theta2[f] - 360.0
        theta3[f]  = degrees(theta3[f]%(2.0*PI))%360.0 
        if theta3[f] >= 180.0:
                theta3[f] = theta3[f] - 360.0
        theta4[f]  = degrees(theta4[f]%(2.0*PI))%360.0 
        if theta4[f] >= 180.0:
                theta4[f] = theta4[f] - 360.0
        theta5[f]  = degrees(theta5[f]%(2.0*PI))%360.0 
        if theta5[f] >= 180.0:
                theta5[f] = theta5[f] - 360.0
        theta6[f]  = degrees(theta6[f]%(2.0*PI))%360.0  
        if theta6[f] >= 180.0:
                theta6[f] = theta6[f] - 360.0

        delta_theta1[f] = theta1[f] - theta2[f]
        delta_theta2[f] = theta1[f] - theta3[f] 
        delta_theta3[f] = theta1[f] - theta4[f] 
        delta_theta4[f] = theta1[f] - theta5[f] 
        delta_theta5[f] = theta1[f] - theta6[f] 

    theta = [theta1,theta2,theta3,theta4,theta5,theta6]

    delta_peri = (per1[0:range_max] - per2[0:range_max] )

    for f in range(0,range_max):

        if delta_peri[f] >= 360.0:
            delta_peri[f] = delta_peri[f]- 360.0
        if delta_peri[f] >= 180.0:
            delta_peri[f] = delta_peri[f]- 360.0
        if delta_peri[f] <= -180.0:
            delta_peri[f] = delta_peri[f]+ 360.0

    len_a = len(a_1)
 

    u1 = G*((1.9891e30)*(Stellar_mass+par[7*0 ]) )
    u2 = G*((1.9891e30)*(Stellar_mass+par[7*0]+par[7*1]) )

    P1 = (((a_1**3.0)*(4*(pi**2)))/u1)**0.5 
    P2 = (((a_2**3.0)*(4*(pi**2)))/u2)**0.5 

    Prat = P2[0:len_a]/P1[0:len_a]
    
    PPrat = mean(Prat[len_a-100:len_a])                     

    if PPrat < 5.98:
        os.chdir("../")  
        os.system("rm -r %s"%rand_dir)    
        return -np.inf
    elif PPrat > 6.1:
        os.chdir("../")  
        os.system("rm -r %s"%rand_dir)    
        return -np.inf

    
    del_peri = (max(delta_peri[len_a-100:len_a]) - min(delta_peri[len_a-100:len_a]))/2.0
    
  
    
    theta_1 = (max(theta1[len_a-100:len_a]) - min(theta1[len_a-100:len_a]))/2.0
    theta_2 = (max(theta2[len_a-100:len_a]) - min(theta2[len_a-100:len_a]))/2.0
    theta_3 = (max(theta3[len_a-100:len_a]) - min(theta3[len_a-100:len_a]))/2.0
    theta_4 = (max(theta4[len_a-100:len_a]) - min(theta4[len_a-100:len_a]))/2.0
    theta_5 = (max(theta5[len_a-100:len_a]) - min(theta5[len_a-100:len_a]))/2.0
    theta_6 = (max(theta6[len_a-100:len_a]) - min(theta6[len_a-100:len_a]))/2.0                   
    

    b_prat = 6.01
    
    b_t1 = 105.3
    b_t2 = 86.0
    b_t3 = 66.8
    b_t4 = 47.6
    b_t5 = 28.9
    b_t6 = 37.3

    d_prat = 0.05
            
    d_t1 = 25.0
    d_t2 = 25.0
    d_t3 = 25.0
    d_t4 = 25.0
    d_t5 = 25.0
    d_t6 = 25.0    
    
    t_o = np.array([ PPrat,  theta_1,theta_2,theta_3,theta_4,theta_5,theta_6])
    t_c = np.array([ b_prat,  b_t1,b_t2,b_t3,b_t4,b_t5,b_t6]        )
    t_s = np.array([ d_prat,  d_t1,d_t2,d_t3,d_t4,d_t5,d_t6]        )
    
    loglik = -0.5*(np.sum(np.divide((t_o - t_c)**2,(d_t1**2))) )
    

#    PPP = Prat[0:range_max]

#        print(loglik,PPrat, theta_1,theta_2,theta_3,theta_4,theta_5,theta_6)
     


#        os.system("python ../plot_migration.py")
#    if PPrat < 5.98:
#                return loglik - 100.0


    os.chdir("../")  
    os.system("rm -r %s"%rand_dir)        
    
    
    f=open("loglik","a")
    f.write("%s %s  %s %s %s %s %s %s %s %s %s %s %s %s %s\n" %(loglik, PPrat, theta_1,theta_2,theta_3,theta_4,theta_5,theta_6, par[7*npl+3],   
    str(par[7*0 +2]),str(par[7*0 +3]),str(par[7*0 +4]),        str(par[7*1 +2]),str(par[7*1 +3]),str(par[7*1 +4])
    ))

    f.close() 

    return loglik
 




############################# Here we start ####################################

npl = 2        # number of planets

mod = "kep" 
#mod = "dyn"
epoch = 2450000.0
output_files = 0

##################### innitial guess of parametrs to vary ######################



# this is really not a flexable format! Will be fixed later on... 

#pl = [K m/s, P days, ecc, peri [deg], Mean_an [deg]]


#pl1 = [563.80278, 256.21505,  0.42615, 161.41258, 182.67848,  90.0, 0.0] # planet 1
#pl2 = [41.10468, 1505.96361, 0.0      , 92.47367, 45.76581, 90.0, 0.0] # planet 2

#  m     a     e     i     omega capom M     rpl   rh
#  3.0991E-01  1.0000E+00  5.0000E-02  0.0000E+00  0.0000E+00  0.0000E+00  1.5708E+00  1.0364E-03  1.3780E-01
#  3.4465E-01  3.5400E+00  5.0000E-02  0.0000E+00  0.0000E+00  0.0000E+00  4.7124E+00  1.0738E-03  5.0540E-01
 
pl1 = [7.85e-3, 1,    0.001, 0.0, 0.0,  0.0, 0.0] # planet 1
pl2 = [8.73e-3, 3.54, 0.12, 0.0, 180.0,  0.0, 0.0] # planet 2

pl1_str = [r'm$_b$ [m/s]',r'a$_b$ [day]' ,r'e$_b$',r'$\omega_p$ [deg]',r'M$_b$ [deg]', r'i_b [deg]',r'\Omega_b [deg]']
pl2_str = [r'm$_c$ [m/s]',r'a$_c$ [day]' ,r'e$_c$',r'$\omega_c$ [deg]',r'M$_c$ [deg]', r'i_c [deg]',r'\Omega_c [deg]'] 
 
use_pl1 = [0, 0, 1, 1, 1, 0, 0] # planet 1
use_pl2 = [0, 0, 1, 1, 1, 0, 0] # planet 2
 
lim_pl1 = [[0.0, 1.0],[0.0, 10.0],[0.0, 0.20],[0, 360.0],[0, 360.0],[3.0, 177.0],[0.0, 359.9]] # planet 1
lim_pl2 = [[0.0, 1.0],[0.0, 10.0],[0.0, 0.20],[0, 360.0],[0, 360.0],[3.0, 177.0],[0.0, 359.9]] # planet 2
 


de = [0, 0 ]
use_de = [0,0]
lim_de = [[0,200],[0,200]]
de_str = ['de1', 'de2' ]


da = [0, -4.8]
use_da = [0,1]
lim_da = [[0,200],[-5.5,-3.9]]
da_str = ['da1', 'da2' ]

 

stmass = [1.00 ] # stellar mass
use_stmass= [0]
lim_stmass= [[0.1,2.5]]
stmass_str = [r'$M_\odot$']


#---------- parameters  -----#
par    = np.concatenate((pl1,pl2,de,da,stmass), axis=0)  
#---------- which parameters to be optimized  -----#
flag   = np.concatenate((use_pl1,use_pl2, use_de,use_da, use_stmass), axis=0)
#------------------------- boundaries for the 'TNC' method --------------------#
bounds = np.concatenate((lim_pl1,lim_pl2,lim_de,lim_da, lim_stmass), axis=0)
#------------------------- Parameter names --------------------#
el_str = np.concatenate((pl1_str, pl2_str,de_str,da_str, stmass_str), axis=0)
#------------------------- Sort out only those that will be optimized --------------------#
p,f,b,e = concentrate_par(par,flag,bounds,el_str)

 

p = []  #'p" are the fitted parameters
b = []
e = []

for j in range(len(par)):
    if flag[j] > 0:
        p.append(par[j])
        b.append(bounds[j])
        e.append(el_str[j])

####### find the -LogLik "minimum" using the "truncated Newton" method ######### 
nll = lambda *args: -migration(*args)

minimzers = ['Nelder-Mead','Powell','CG','BFGS','Newton-CG','L-BFGS-B',
'TNC','COBYLA','SLSQP','dogleg','trust-ncg']

for k in range(0): # run at least 3 times the minimizer

    result = op.minimize(nll, p, args=(par,f,  
npl), method=minimzers[6],bounds=b, options={'xtol': 1e-6, 
'disp': True })

    p = result["x"]
    print("Best fit par.:", p)
#----------------- one more time using the Simplex method ---------------------#

for k in range(0): # run at least 3 times the minimizer
    result = op.minimize(nll, result["x"], args=(par,f,  
npl), method=minimzers[0], options={'xtol': 1e-6, 'disp': True,
 'maxiter':30000, 'maxfev':30000 })


################################################################################ 

best_fit_par = p# result["x"]
print("Best fit par.:", best_fit_par)



####################### Start the mcmc using "emcee"   ######################### 

#------------------------- flat prior for now!!! ------------------------------#

b = np.array(b)

 
def lnprior(p): 
    for j in range(len(p)):
 
 
        ######## if something is away from the borders - reject #####
        if p[j] <= b[j,0] or p[j] >= b[j,1]:
            return -np.inf
        ######## normal distribution for the eccentricity with
        ######## scipy.stats.norm(mu,std_dev).pdf(e[j]) 
        ######## flat prior for the rest! 
        
       # if j == 3:
            #prob_ecc1 =  np.log(float(scipy.stats.norm(0.001,0.02).pdf(p[j])) )
       #     prob_ecc1 = np.log(np.exp(-0.5 * (p[j]/SIG_E)**2.0))

       # if j == 8:

        #    prob_ecc2 = np.log(np.exp(-0.5 * (p[j]/SIG_E)**2.0))

    return 0.0 # + prob_ecc1 + prob_ecc2
  

####################################################### 

def lnprob(p, stmass, npl):
    lp = lnprior(p)
    if not np.isfinite(lp):
        return -np.inf
    return lp + migration(p, par,f, npl)

 

################## prior TESTS ########################

def prior_transform(p): 
 
        u_trans = np.zeros(len(p)) 
        for jj in range(len(p)): 

                #            if priors[0][j,2] == True:
                #                u_trans[j] = trans_norm(p[j],bb[j][0],bb[j][1])
                #            elif priors[1][j,2] == True:
                #                u_trans[j] = trans_loguni(p[j],bb[j][0],bb[j][1]) 
                #            else:
                u_trans[jj] = trans_uni(p[jj],b[jj][0],b[jj][1])
        return u_trans   


def trans_norm(p ,mu,sig):
    return stats.norm.ppf(p,loc=mu,scale=sig)

def trans_uni(p,a,b):
    return a + (b-a)*p

def trans_loguni(p,a,b):
    return np.exp(np.log(a) + p*(np.log(b)-np.log(a)))    
    

################## partial_func ########################

partial_func = FunctionWrapper(migration, (par,f, npl) )




print(p)
#priors = [pr_nr,jeff_nr]

ndim, nwalkers = len(p), len(p)*5
 
print_progress = True
Dynamic_nest = True
threads = 8
dynesty_samp = 'rwalk'
stop_crit = 0.1
ns_bound = 'multi'
ns_pfrac = 1.0
ns_use_stop = True

ns_maxiter = 100000
ns_maxcall = 500

nest_weighted = False
fileoutput = True


################## NS sampler #########################

from pathos.pools import ProcessPool as Pool
 
 

pool = Pool(ncpus=threads)#, context=ctx)
#thread = Pool(threads)#, context=ctx)



if Dynamic_nest == False:
    print("'Static' Nest. Samp. is running, please wait...")

    if threads > 1:

        sampler = dynesty.NestedSampler(partial_func, prior_transform, ndim, nlive=nwalkers, pool = pool,
                                            queue_size=threads, sample = dynesty_samp, bound = ns_bound)
        sampler.run_nested(print_progress=print_progress,dlogz=stop_crit, 
        maxiter = ns_maxiter, maxcall = ns_maxcall ) #dlogz=stop_crit,
        #pool.close()
        #pool.join()
        #pool.clear()

    else:
         sampler = dynesty.NestedSampler(partial_func, prior_transform, ndim, nlive=nwalkers, sample = dynesty_samp, bound = ns_bound)
         sampler.run_nested(print_progress=print_progress,dlogz=stop_crit, 
         maxiter = ns_maxiter, maxcall = ns_maxcall )


    #obj.dyn_res = sampler.results
    sampler.results.summary()

else:
    print("'Dynamic' Nest. Samp. is running, please wait...")

    if threads > 1:

   
        sampler = dynesty.DynamicNestedSampler(partial_func, prior_transform, ndim, pool = pool,
                                               queue_size=threads, sample = dynesty_samp, bound = ns_bound) # nlive=nwalkers,

        sampler.run_nested(print_progress=print_progress,dlogz_init=stop_crit,nlive_init=nwalkers, 
        maxiter = ns_maxiter, maxcall = ns_maxcall,use_stop = ns_use_stop, wt_kwargs={'pfrac': ns_pfrac})   #nlive_batch=1
        #pool.close()
        #pool.join()
        #pool.clear()
           
        

    else:
        sampler = dynesty.DynamicNestedSampler(partial_func, prior_transform, ndim, sample = dynesty_samp, bound = ns_bound)
        sampler.run_nested(print_progress=print_progress,dlogz_init=stop_crit,nlive_init=nwalkers,  
        maxiter = ns_maxiter, maxcall = ns_maxcall,use_stop = ns_use_stop, wt_kwargs={'pfrac': ns_pfrac} ) 

    # just in case
    pool.close()
    pool.join()
    pool.clear()

   # obj.dyn_res = sampler.results

    res = ("niter: {:d}\n"
           "ncall: {:d}\n"
           "eff(%): {:6.3f}\n"
           "logz: {:6.3f} +/- {:6.3f}".format(sampler.results.niter, sum(sampler.results.ncall),
                   sampler.results.eff, sampler.results.logz[-1], sampler.results.logzerr[-1]))

    print('Summary\n=======\n'+res)



# print("--- %s seconds ---" % (time.time() - start_time))



ln = np.hstack(sampler.results.logl)



if nest_weighted == True:
    weighted = np.exp(sampler.results.logwt - sampler.results.logz[-1])
    samples  =  dill.copy(dynesty.utils.resample_equal(sampler.results.samples, weighted))   
else:
    samples  =  dill.copy(sampler.results.samples)   


# print("--- %s seconds ---" % (time.time() - start_time))  




if (fileoutput):
# start_time = time.time()   
# print("Please wait... writing the ascii file")  

    outfile = open(str("nest_sampl"), 'w') # file to save samples
    for j in range(len(samples)):
        outfile.write("%s  " %(ln[j]))
        for z in range(len(p)):
            outfile.write("%s  " %(samples[j,z]))
        outfile.write("\n")
    outfile.close()        
# print("--- Done for ---")           
#   print("--- %s seconds ---" % (time.time() - start_time))  
 
#start_time = time.time()    

############### find uncertainties form the result distribution#################

#----------------------------------- labels  ----------------------------------#
 
cornerplot_opt = {'bins': 25,
  'color': 'k',
  'reverse': False,
  'upper': True,
  'quantiles': 68.3,
  'levels': (0.6827, 0.9545, 0.9973),
  'smooth': 1.0,
  'smooth1d': 1.0,
  'plot_contours': True,
  'show_titles': False,
  'dpi': 300,
  'pad': 15,
  'labelpad': 0.09,
  'truth_color': 'r',
  'title_kwargs': {'fontsize': 12},
  'scale_hist': True,
  'fill_contours': True,
  'no_fill_contours': False,
  'plot_datapoints': False,
  'stab_color': 'r',
  'stab_threshold': 1.0}
 
 
level = (100.0-68.3)/2.0
#quantiles = None
#title_quantiles = None    
level_q = (100.0-cornerplot_opt["quantiles"])/2.0    

quantiles = [level_q/100.0, 0.5, 1.0-level_q/100.0]
title_quantiles = [level_q/100.0, 0.5, 1.0-level_q/100.0]  



best_fit_par_2 = []
print("Best fit par. and their 1 sigma errors" )	
for i in range(len(best_fit_par)):
    ci = np.percentile(samples[:,i], [level, 100.0-level])
    print(e[i],'=', best_fit_par[i], "- %s"%(best_fit_par[i]-ci[0]), "+ %s"%(ci[1]  - best_fit_par[i] ))

print("   ")
print("   " 	)
print("Means and their 1 sigma errors" )	
for i in range(len(best_fit_par)):
    ci = np.percentile(samples[:,i], [level, 100.0-level])
    print(e[i],'=', np.mean(samples[:,i]), "- %s"%(np.mean(samples[:,i])-ci[0]), "+ %s"%(ci[1]  - np.mean(samples[:,i]) ))

    best_fit_par_2.append(np.median(samples[:,i]))


#fig = corner.corner(samples, labels=e, truths=best_fit_par, dpi = 300 )
#fig.savefig("samples_%s.png"%mod)




 
#range=ranged, 
#fig = corner.corner(samples,bins=25, color="k", reverse=True, upper= True, labels=e, quantiles=[0.1585, 0.8415],levels=(0.6827, 0.9545,0.9973), smooth=1.0, smooth1d=1.0, plot_contours= True, show_titles=True, truths=best_fit_par, dpi = 300, pad=15, labelpad = 50 ,truth_color ='r', title_kwargs={"fontsize": 12}, scale_hist=True,  no_fill_contours=True, plot_datapoints=True)

#fig = corner.corner(samples,bins=25, color="k", reverse=False, upper= True, labels=e, 
#                    quantiles=[0.1585, 0.8415],levels=(0.6827, 0.9545,0.9973), smooth=1.0, 
#                    smooth1d=1.0, plot_contours= True, show_titles=True, truths=best_fit_par_2, dpi = 300, 
#                    pad=15, labelpad = 0 ,truth_color ='r', title_kwargs={"fontsize": 12}, scale_hist=True,  
#                    no_fill_contours=True, plot_datapoints=True)

fig = corner.corner(
    samples,
    bins=cornerplot_opt["bins"],
    color=cornerplot_opt["color"],
    reverse=cornerplot_opt["reverse"],
    upper=cornerplot_opt["upper"],
    labels=e,
    quantiles=quantiles,
    title_quantiles=title_quantiles,
    levels=(0.6827, 0.9545, 0.9973),
    smooth=cornerplot_opt["smooth"],
    smooth1d=cornerplot_opt["smooth1d"],
    plot_contours=cornerplot_opt["plot_contours"],
    show_titles=cornerplot_opt["show_titles"],
    truths=best_fit_par,
    dpi=cornerplot_opt["dpi"],
    pad=cornerplot_opt["pad"],
    labelpad=cornerplot_opt["labelpad"],
    truth_color=cornerplot_opt["truth_color"],
    title_kwargs={"fontsize": 12},
    scale_hist=cornerplot_opt["scale_hist"],
    no_fill_contours=cornerplot_opt["no_fill_contours"],
    fill_contours=cornerplot_opt["fill_contours"],
    plot_datapoints=cornerplot_opt["plot_datapoints"],
    contour_kwargs={'colors': cornerplot_opt["stab_color"]},
    hist_kwargs={'color': cornerplot_opt["stab_color"]},
    data_kwargs={'zorder': 10, 'color': cornerplot_opt["stab_color"]}
)


#fig = corner.corner(samples, labels=el_str, truths=best_fit_par, dpi = 300, pad=15, labelpad = 50 )
fig.savefig("samples_%s.pdf"%mod)




os.system("sort nest_sampl | uniq > sorted_nest_sampl &")

print("Done")
















