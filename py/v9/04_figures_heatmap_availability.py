# 04_figures_heatmap_availability.py
# v9 rewrite (2026-11-08): Figure 2 (type I error heat map), Figure S1 (usable share heat map) and Figure 3 (availability and power).
# Cells extracted verbatim, in execution order, from the session that produced them (see py/v9/provenance/);
# only file paths were changed (artifact store -> tables/v9/). Run from the repository root.
# Requires the helpers apply_figure_style() and panel_letter() of the figure-style library used in that session (not part of this repository)
# and the objects defined by 02_results_by_mechanism.py (run it first in the same Python session).

# ---- cell 18: Drawing heat-map figures A and B
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, Normalize
from matplotlib.patches import Rectangle
apply_figure_style(sizes=(7,6,6))
cols=[(t,n) for t in Tb for n in Nb]
nlab={'N=10':'10','N=15-25':'15–25','N=35-65':'35–65','N=75-145':'75–145'}
mtitle={'GenBinaryY':'GenBinaryY (latent + transition)','GLMM-native':'GLMM-native (random intercept)','GLMM-AR1':'GLMM-AR1 (random intercept + AR(1))','Copula-AR1':'Copula-AR1 (serial correlation only)'}
R7=reg[reg.strategy.isin(main7)].copy(); R7['strategy']=R7.strategy.astype(str)
def mat(m,var):
    g=R7[R7.dgp==m].set_index(['strategy','Tband','Nband'])
    return np.array([[g.loc[(s,t,n),var] for (t,n) in cols] for s in main7])
def lum(rgba): r,g_,b=rgba[:3]; return 0.2126*r+0.7152*g_+0.0722*b
def heat(var,fname,cmap,norm,fmt,cbar_label,outline=True):
    fig,axs=plt.subplots(2,2,figsize=(7,4.9),gridspec_kw=dict(hspace=0.42,wspace=0.04,left=0.17,right=0.985,top=0.9,bottom=0.12))
    for i,(ax,m) in enumerate(zip(axs.flat,mechs)):
        V=mat(m,var); A=mat(m,'n_above'); NA=mat(m,'n_assessed')
        Vc=np.clip(V,norm.vmin,norm.vmax)
        im=ax.imshow(np.ma.masked_invalid(Vc),cmap=cmap,norm=norm,aspect='auto')
        ax.set_facecolor('0.85')
        for r in range(V.shape[0]):
            for c in range(V.shape[1]):
                v=V[r,c]
                if np.isnan(v): ax.text(c,r,'–',ha='center',va='center',fontsize=6,color='0.3'); continue
                col=cmap(norm(np.clip(v,norm.vmin,norm.vmax)))
                ax.text(c,r,fmt(v),ha='center',va='center',fontsize=6,color='white' if lum(col)<0.45 else 'black',
                        fontweight='bold' if (outline and A[r,c]>0) else 'normal')
                if outline and A[r,c]>0:
                    ax.add_patch(Rectangle((c-0.5+0.04,r-0.5+0.04),0.92,0.92,fill=False,lw=1.3,ec='black'))
        for x in [3.5,7.5]: ax.axvline(x,color='white',lw=2.5)
        ax.set_xticks(range(12)); ax.set_xticklabels([nlab[n] for t,n in cols],rotation=90,fontsize=6)
        ax.set_yticks(range(7)); ax.set_yticklabels([labels[s] for s in main7] if i%2==0 else [],fontsize=6)
        ax.tick_params(length=0)
        for s in ax.spines.values(): s.set_visible(False)
        for j,t in enumerate(Tb):
            ax.text(1.5+4*j,-0.75,t,ha='center',va='bottom',fontsize=6)
        ax.set_title(mtitle[m],loc='left',fontsize=7,pad=11)
        if i>=2: ax.set_xlabel('Number of clusters N',fontsize=6,labelpad=2)
        panel_letter(ax,'abcd'[i])
    cax=fig.add_axes([0.35,0.035,0.45,0.018])
    cb=fig.colorbar(mpl.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax,orientation='horizontal')
    cb.set_label(cbar_label,fontsize=6); cb.ax.tick_params(labelsize=6)
    fig.savefig(f'tables/v9/{fname}.pdf'); fig.savefig(f'tables/v9/{fname}.png',dpi=300)
    return fig
cmapA=plt.get_cmap('RdBu_r')
normA=TwoSlopeNorm(vmin=0,vcenter=5,vmax=15)
figA=heat('mean_rate','figA_typeI_heatmap',cmapA,normA,lambda v:f'{v:.1f}','Mean type I error (%) over assessed design points; nominal 5% = white; clipped at 0 and 15. Boxed, bold: ≥1 point clearly above 7.5%')
cmapB=plt.get_cmap('viridis')
figB=heat('mean_pu','figB_usable_heatmap',cmapB,Normalize(0,100),lambda v:f'{v:.0f}','Usable share under H0 (%) — mean over design points',outline=False)
print('ok')

# ---- cell 19: Fixing heat-map layout overlaps
def heat(var,fname,cmap,norm,fmt,cbar_label,outline=True,ticks=None):
    fig,axs=plt.subplots(2,2,figsize=(7,5.4),gridspec_kw=dict(hspace=0.5,wspace=0.04,left=0.17,right=0.985,top=0.9,bottom=0.155))
    for i,(ax,m) in enumerate(zip(axs.flat,mechs)):
        V=mat(m,var); A=mat(m,'n_above')
        ax.imshow(np.ma.masked_invalid(np.clip(V,norm.vmin,norm.vmax)),cmap=cmap,norm=norm,aspect='auto')
        ax.set_facecolor('0.85')
        for r in range(V.shape[0]):
            for c in range(V.shape[1]):
                v=V[r,c]
                if np.isnan(v): ax.text(c,r,'n.a.',ha='center',va='center',fontsize=6,color='0.25'); continue
                col=cmap(norm(np.clip(v,norm.vmin,norm.vmax)))
                flag=outline and A[r,c]>0
                ax.text(c,r,fmt(v),ha='center',va='center',fontsize=6,color='white' if lum(col)<0.45 else 'black',fontweight='bold' if flag else 'normal')
                if flag: ax.add_patch(Rectangle((c-0.46,r-0.46),0.92,0.92,fill=False,lw=1.3,ec='black'))
        for x in [3.5,7.5]: ax.axvline(x,color='white',lw=2.5)
        ax.set_xticks(range(12)); ax.set_xticklabels([nlab[n] for t,n in cols],rotation=90,fontsize=6)
        ax.set_yticks(range(7)); ax.set_yticklabels([labels[s] for s in main7] if i%2==0 else [],fontsize=6)
        ax.tick_params(length=0)
        for s in ax.spines.values(): s.set_visible(False)
        for j,t in enumerate(Tb): ax.text(1.5+4*j,-0.7,t,ha='center',va='bottom',fontsize=6)
        ax.text(0,1.2,'abcd'[i],transform=ax.transAxes,fontsize=8,fontweight='bold',ha='left',va='bottom')
        ax.text(0.045,1.2,mtitle[m],transform=ax.transAxes,fontsize=7,ha='left',va='bottom')
        if i>=2: ax.set_xlabel('Number of clusters N',fontsize=6,labelpad=2)
    cax=fig.add_axes([0.3,0.045,0.5,0.016])
    cb=fig.colorbar(mpl.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax,orientation='horizontal')
    if ticks is not None: cb.set_ticks(ticks)
    cb.set_label(cbar_label,fontsize=6); cb.ax.tick_params(labelsize=6)
    fig.savefig(f'tables/v9/{fname}.pdf'); fig.savefig(f'tables/v9/{fname}.png',dpi=300)
    r_=fig.canvas.get_renderer()
    T=[(t,t.get_window_extent(r_)) for t in fig.findobj(mpl.text.Text) if t.get_text().strip() and t.get_visible()]
    ov=[(a.get_text(),b.get_text()) for k,(a,ba) in enumerate(T) for b,bb in T[k+1:] if ba.overlaps(bb)]
    print(fname,'overlaps:',ov[:10])
    return fig
figA=heat('mean_rate','figA_typeI_heatmap',cmapA,normA,lambda v:f'{v:.1f}','Mean type I error (%) over assessed points (5% = white; clipped 0–15). Boxed bold cell: ≥1 design point clearly above 7.5%',ticks=[0,2.5,5,7.5,10,15])
figB=heat('mean_pu','figB_usable_heatmap',cmapB,Normalize(0,100),lambda v:f'{v:.0f}','Usable share of datasets under H0 (%), mean over design points',outline=False,ticks=[0,25,50,75,100])

# ---- cell 20: Drawing figure C availability and power
M7=M[M.strategy.isin(main7)].copy(); M7['strategy']=M7.strategy.astype(str)
byT=M7.groupby(['strategy','ntime']).agg(mean_pu=('pu','mean'),mean_upow=('upow','mean'),n_points=('pu','size')).reset_index()
byT.round(2).to_csv('tables/v9/figC_data_byT.csv',index=False)
print(sorted(M7.ntime.unique()), byT.n_points.unique())
oi=['#000000','#E69F00','#56B4E9','#009E73','#D55E00','#0072B2','#CC79A7']
mk=['o','s','^','v','D','P','X']; ls=['-','-','-','--','-','-',':']
fig,axs=plt.subplots(1,2,figsize=(7,3.0),gridspec_kw=dict(wspace=0.28,left=0.08,right=0.99,top=0.86,bottom=0.33))
for ax,var,yl,tt in [(axs[0],'mean_pu','Usable share under H0 (%)','Usable share rises with T for every analysis'),
                     (axs[1],'mean_upow','Power per attempted dataset (%)','Power per attempted dataset, H1')]:
    for s,c,m_,l in zip(main7,oi,mk,ls):
        g=byT[byT.strategy==s]
        ax.plot(g.ntime,g[var],color=c,marker=m_,ls=l,ms=3.5,lw=1.2,label=labels[s])
    ax.set_xlabel('Number of time points T'); ax.set_ylabel(yl); ax.set_xticks(sorted(M7.ntime.unique()))
    ax.margins(0.04)
axs[0].set_ylim(0,102)
panel_letter(axs[0],'a'); panel_letter(axs[1],'b')
axs[0].set_title('Usable share by T',loc='left'); axs[1].set_title('Power by T',loc='left')
h,l_=axs[0].get_legend_handles_labels()
fig.legend(h,l_,loc='lower center',ncol=4,frameon=False,fontsize=6,bbox_to_anchor=(0.53,0.0))
fig.savefig('tables/v9/figC_availability_power.pdf'); fig.savefig('tables/v9/figC_availability_power.png',dpi=300)
print(byT.pivot(index='strategy',columns='ntime',values='mean_pu').round(1).loc[main7].to_string())
print(byT.pivot(index='strategy',columns='ntime',values='mean_upow').round(1).loc[main7].to_string())

# ---- cell 21: Polishing figure C and saving
axs[1].set_ylim(0,60); fig.subplots_adjust(bottom=0.27)
axs[1].set_title('Power by T (uncorrected GEE and GLMM Wald are anti-conservative)',loc='left',fontsize=6)
fig.savefig('tables/v9/figC_availability_power.pdf'); fig.savefig('tables/v9/figC_availability_power.png',dpi=300); print('ok')

