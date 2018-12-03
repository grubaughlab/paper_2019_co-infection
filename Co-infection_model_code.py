# -*- coding: utf-8 -*-
"""
Created on Tue Aug 14 11:16:09 2018

@author: sean
"""

import scipy.integrate as spi
import numpy as np
import matplotlib.pyplot as plt
from ternary.helpers import simplex_iterator
import ternary
import matplotlib.image as mpimg

plt.ioff()

#parameter values:
m = 1.0
a = 0.5
b1 = b2 = b12 = b = 0.3
c1 = c2 = c12 = 0.5
alpha1H = alpha2H = alpha1M = alpha2M = 0.0
r = 5.0
g = 1.0/8
ND = 1.75*365.25
TS = 1.0
virus2introduction = 30.0

def diff_eqs(INP,t):  
    '''The main set of equations'''
    [SH, I1H, I2H, I12CH, I12SH, I1_2H, I2_1H, R1H, R2H, R12H,\
    SM, I1M, I2M, I12M, cumulative_I12CH, cumulative_I12SH] = INP
    lambda1H = m*a*b1*(I1M+I12M*p1H)
    lambda2H = m*a*b2*(I2M+I12M*p2H)
    lambda12H = m*a*b12*I12M*p12H
    
    lambda1M = a*c1*(I1H+I2_1H+(I12CH+I12SH)*p1M)
    lambda2M = a*c2*(I2H+I1_2H+(I12CH+I12SH)*p2M)
    lambda12M = a*c12*(I12CH+I12SH)*p12M
        
    dSH = -(lambda1H+lambda2H+lambda12H)*SH
    dI1H = lambda1H*SH - ((lambda2H + lambda12H)*(1-alpha1H) + 1./r)*I1H
    dI2H = lambda2H*SH - ((lambda1H + lambda12H)*(1-alpha2H) + 1./r)*I2H
    dI12CH = lambda12H*SH + lambda12H*((1-alpha2H)*I2H + (1-alpha1H)*I1H) - (1./r)*I12CH
    dI12SH = lambda2H*(1-alpha1H)*I1H + lambda1H*(1-alpha2H)*I2H - (1./r)*I12SH
    
    dI1_2H = (lambda12H+lambda2H)*(1-alpha1H)*R1H - (1./r)*I1_2H
    dI2_1H = (lambda12H+lambda1H)*(1-alpha2H)*R2H - (1./r)*I2_1H
    dR1H = (1./r)*I1H - (lambda12H+lambda2H)*(1-alpha1H)*R1H
    dR2H = (1./r)*I2H - (lambda12H+lambda1H)*(1-alpha2H)*R2H
    dR12H = (1./r)*(I2_1H + I1_2H + I12SH + I12CH)
    
    dSM = g*(I1M+I2M+I12M) - (lambda1M+lambda2M+lambda12M)*SM
    dI1M = lambda1M*SM - ((lambda2M+lambda12M)*(1-alpha1M) + g)*I1M
    dI2M = lambda2M*SM - ((lambda1M+lambda12M)*(1-alpha2M) + g)*I2M
    dI12M = lambda12M*SM + (lambda2M+lambda12M)*(1-alpha1M)*I1M + (lambda1M+lambda12M)*(1-alpha2M)*I2M - g*I12M
    
    dcumulative_I12CH = lambda12H*SH + lambda12H*((1-alpha2H)*I2H + (1-alpha1H)*I1H)
    dcumulative_I12SH = lambda2H*(1-alpha1H)*I1H + lambda1H*(1-alpha2H)*I2H
    
    Y = np.array([dSH, dI1H, dI2H, dI12CH, dI12SH, dI1_2H, dI2_1H, dR1H, dR2H, dR12H, \
    dSM, dI1M, dI2M, dI12M, dcumulative_I12CH, dcumulative_I12SH])

    return Y   # For odeint





'''overall dynamics with and without cotransmission - Fixed b'''
#one inital mosquito with virus X (if population a million)
INPUT = np.zeros(16)
INPUT[1] = 1e-6
INPUT[0] = 1-sum(INPUT[1:10])
INPUT[10] = 1-sum(INPUT[11:14])

#parameters for human-> mosquito cotransmission - 60%
p1M = p2M = 0.2
p12M = 1-p1M-p2M
#parameters for mosquito -> human cotransmission - 50% initially, varied later.
p1H = p2H = 0.25
p12H = 1-p1H-p2H

#Times
t_start = 0.0; t_end = ND; t_inc = TS
#this is the time that just virus X is circulating:
t_range1 = np.arange(t_start, virus2introduction, t_inc) 
#this is the time that both viruses are circulating:
t_range2 = np.arange(virus2introduction, t_end+t_inc, t_inc)
t_range = np.concatenate((t_range1, t_range2))

#Set up patterns and labels for figures
colors=['#1A9099', '#A6BDDB', '#E6605F']
hatch_colors = ['#E6605F','#F0B0AF']
hatching = ['+', '.']
linestyles = ['-', '--', '-.']
labels = ["Virus X only", "Virus Y only", "Co-infection"]
hatch_labels = ["Co-transmission", "Sequential infection"]
b1 = b2 = b12 = b
    
#50% co-transmission
RES1 = spi.odeint(diff_eqs,INPUT,t_range1)
#Import a second infection
INPUT2 = RES1[-1] 
INPUT2[2] = 1e-6 
RES2 = spi.odeint(diff_eqs, INPUT2, t_range2)
RES = np.concatenate((RES1, RES2))
[SH, I1H, I2H, I12CH, I12SH, I1_2H, I2_1H, R1H, R2H, R12H,\
    SM, I1M, I2M, I12M, cumulative_I12CH, cumulative_I12SH] = RES.T
#virus X, virus Y and co-infection, stacked, converted to population size:
y = np.row_stack(((I1H+I2_1H), (I2H+I1_2H),(I12CH+I12SH)))
O2 = np.int_(y*1e6+0.5)
#co-transmission and sequential transmission stacked, converted to population size:
yC = np.row_stack((I12CH, I12SH))
C2 = np.int_(yC*1e6+0.5)

#No cotransmission    
p12H = 0
p1H=p2H=0.5
RES1 = spi.odeint(diff_eqs,INPUT,t_range1)
#Import a second infection
INPUT2 = RES1[-1]
INPUT2[2] = 1e-6
RES2 = spi.odeint(diff_eqs, INPUT2, t_range2)
RES = np.concatenate((RES1, RES2))
[SHB, I1HB, I2HB, I12CHB, I12SHB, I1_2HB, I2_1HB, R1HB, R2HB, R12HB,\
    SMB, I1MB, I2MB, I12MB, cumulative_I12CHB, cumulative_I12SHB] = RES.T
#virus X, virus Y and co-infection, stacked, converted to population size:
y = np.row_stack(((I1HB+I2_1HB), (I2HB+I1_2HB),(I12CHB+I12SHB)))
O1 = np.int_(y*1e6+0.5)
#co-transmission and sequential transmission stacked, converted to population size:
yC = np.row_stack((I12CHB, I12SHB))
C1 = np.int_(yC*1e6+0.5)

#100% cotransmission
p1H = p2H = 0.
p12H = 1-p1H-p2H
RES1 = spi.odeint(diff_eqs,INPUT,t_range1)
#Import a second infection
INPUT2 = RES1[-1]
INPUT2[2] = 1e-6
RES2 = spi.odeint(diff_eqs, INPUT2, t_range2)
RES = np.concatenate((RES1, RES2))
[SH, I1H, I2H, I12CH, I12SH, I1_2H, I2_1H, R1H, R2H, R12H,\
    SM, I1M, I2M, I12M, cumulative_I12CH, cumulative_I12SH] = RES.T
#virus X, virus Y and co-infection, stacked, converted to population size:
y = np.row_stack(((I1H+I2_1H), (I2H+I1_2H),(I12CH+I12SH)))
O3 = np.int_(y*1e6+0.5)
#co-transmission and sequential transmission stacked, converted to population size:
yC = np.row_stack((I12CH, I12SH))
C3 = np.int_(yC*1e6+0.5)

fig, axes = plt.subplots(nrows=3, ncols=4, figsize=(24,16))
fig.patch.set_facecolor('white')

#human mosquito image
im = mpimg.imread('Mosquito_To_Human.png')
axes[0][0].imshow(im, cmap ='gist_gray')
axes[0][0].axis('off')

#Pie charts to indicate co-transmission proportions
s2 = [0.25, 0.25, 0.5]
s1 = [0.5, 0.5, 0.0]
s3 = [0.0,0.0, 1.0]
axes[0][1].pie(s1, colors=colors, startangle=90+int(0.5*s1[-1]*360), radius=0.5)
axes[0][2].pie(s2, colors=colors, startangle=90+int(0.5*s2[-1]*360), radius=0.5)                
axes[0][3].pie(s3, colors=colors, startangle=90+int(0.5*s3[-1]*360), radius=0.5)     
axes[0][1].set_aspect('equal')           
axes[0][2].set_aspect('equal')           
axes[0][3].set_aspect('equal')           
axes[0][1].set_title('No co-transmission', fontsize="xx-large")           
axes[0][2].set_title('50% co-transmission', fontsize="xx-large")                
axes[0][3].set_title('100% co-transmission', fontsize="xx-large")      

#stacked lineplots
axes[1][1].stackplot(np.arange(len(SH)),O1, colors=colors, linestyles=linestyles, labels=labels)
stacks = axes[2][1].stackplot(np.arange(len(SH)),C1, colors=hatch_colors,  labels=hatch_labels)
for stack, hatch in zip(stacks, hatching):
    stack.set_hatch(hatch)
    
axes[1][2].stackplot(np.arange(len(SH)),O2, colors=colors, linestyles=linestyles, labels=labels)
stacks = axes[2][2].stackplot(np.arange(len(SH)),C2, colors=hatch_colors,  labels=hatch_labels)
for stack, hatch in zip(stacks, hatching):
    stack.set_hatch(hatch)

axes[1][3].stackplot(np.arange(len(SH)),O3, colors=colors, linestyles=linestyles, labels=labels)
axes[1][3].axis(sharey=axes[1][1])
stacks = axes[2][3].stackplot(np.arange(len(SH)),C3, colors=hatch_colors,  labels=hatch_labels)
for stack, hatch in zip(stacks, hatching):
    stack.set_hatch(hatch)

#Get and set axis limits on line plots:
xlimits = axes[1][1].get_xlim()
ylimits_overall = (0.0, max(axes[1][1].get_ylim()[1],axes[1][2].get_ylim()[1],axes[1][3].get_ylim()[1]))
ylimits_coinfection = (0.0, max(axes[2][1].get_ylim()[1],axes[2][2].get_ylim()[1],axes[2][3].get_ylim()[1]))
for i in xrange(3):
    for j in xrange(2):
        axes[j+1][i+1].set_xlim(xlimits)
        if j:
            axes[j+1][i+1].set_ylim(ylimits_coinfection)
            axes[j+1][i+1].set_xlabel("Time (days)", fontsize = "x-large")
            axes[j+1][i+1].set_xticklabels(np.delete(np.int_(axes[j+1][i+1].get_xticks()),-1), fontsize="large")
        else:
            axes[j+1][i+1].set_ylim(ylimits_overall)
            axes[j+1][i+1].set_xticklabels([])
        if i:
            axes[j+1][i+1].set_yticklabels([])
        else:
            axes[j+1][i+1].set_ylabel("Number of cases", fontsize = "x-large")
            axes[j+1][i+1].set_yticklabels(np.delete(np.int_(axes[j+1][i+1].get_yticks()),-1), fontsize="large")

#Legend and labels
legend = axes[1][0].legend(*axes[1][1].get_legend_handles_labels(), loc='center left', ncol = 1, fontsize='x-large')
legend.set_title("All cases", prop = {'size':'xx-large'})
legend = axes[2][0].legend(*axes[2][1].get_legend_handles_labels(), loc='center left', ncol = 1, fontsize='x-large')
legend.set_title("Co-infection cases", prop = {'size':'xx-large'})
axes[1][0].axis('off')
axes[2][0].axis('off')

fig.subplots_adjust(wspace=0.01, hspace=0.01)
#fig.savefig("figures/complete_figure.pdf", dpi =300, bbox_inches='tight')
fig.show()





'''barycentric plot of cotransmission probabilities'''
def generate_random_heatmap_data(scale=5):
    d = dict()
    c = dict()
    for (i,j,k) in simplex_iterator(scale):
        global p1H
        global p2H     
        global p12H
        p1H, p2H, p12H = float(i)/scale, float(j)/scale, float(k)/scale
        t_range1 = np.arange(t_start, virus2introduction, t_inc)
        t_range2 = np.arange(virus2introduction, t_end+t_inc, t_inc)
        RES1 = spi.odeint(diff_eqs,INPUT,t_range1)    
        INPUT2 = RES1[-1]
        INPUT2[2] = 1e-6
        RES2 = spi.odeint(diff_eqs, INPUT2, t_range2)
        RES = np.concatenate((RES1, RES2))
        [cumulative_I12CH, cumulative_I12SH] = RES.T[-2:,-1]        
        d[(i,j)] = cumulative_I12CH/(cumulative_I12CH+cumulative_I12SH)
        c[(i,j)] = cumulative_I12CH+cumulative_I12SH
    return d, c

scale = 100
c = generate_random_heatmap_data(scale=scale)[1]
figure, tax = ternary.figure(scale=scale)
tax.heatmap(data=c, scale=scale, cbarlabel='Total proportion co-infected')
tax.boundary()
tax.clear_matplotlib_ticks()
tax.ticks(multiple = 10)
tax.left_axis_label("$p_{12}$ (%)")
tax.right_axis_label("$p_{1}$ (%)")
tax.bottom_axis_label("$p_{2}$ (%)")
tax.ax.axis('off')
#tax.savefig("figures/ternary_total_coinfected_probabilities.pdf", dpi =300)
tax.show()

#revert parameters to initial values
p1H = 0.25
p2H = 0.25
p12H = 0.5





'''line plot of cotransmission probabilities'''
p12Hset = np.arange(0.0,1.01, 0.01)
cumulative_cotransmission = np.zeros(len(p12Hset))
cumulative_sequential = np.zeros(len(p12Hset))
t_range1 = np.arange(t_start, virus2introduction, t_inc)
t_range2 = np.arange(virus2introduction, t_end+t_inc, t_inc)
t_range = np.concatenate((t_range1, t_range2))

i=0
for p12H in p12Hset:
    p1H = p2H = (1-p12H)/2.
    RES1 = spi.odeint(diff_eqs,INPUT,t_range1)    
    INPUT2 = RES1[-1]
    INPUT2[2] = 1e-6
    RES2 = spi.odeint(diff_eqs, INPUT2, t_range2)
    RES = np.concatenate((RES1, RES2))
    cumulative_cotransmission[i] = RES[-1,-2]
    cumulative_sequential[i] = RES[-1,-1]
    i+=1

fig, (ax1, ax2) = plt.subplots(1,2)  
ax1.plot(p12Hset, cumulative_cotransmission+cumulative_sequential) 
ax2.plot(p12Hset, cumulative_cotransmission/(cumulative_cotransmission+cumulative_sequential)) 
ax1.set_xlabel('probability of co-transmission')
ax2.set_xlabel('probability of co-transmission')
ax1.set_ylabel('prevalence of coinfection')
ax2.set_ylabel('proportion due to cotransmission')    
ax1.set_aspect((ax1.get_xlim()[1]-ax1.get_xlim()[0])/(ax1.get_ylim()[1]-ax1.get_ylim()[0]))
ax2.set_aspect((ax2.get_xlim()[1]-ax2.get_xlim()[0])/(ax2.get_ylim()[1]-ax2.get_ylim()[0])) 

fig.tight_layout()
#fig.savefig("figures/line_plot_co_transmission_probability.pdf", dpi =300)
fig.show()