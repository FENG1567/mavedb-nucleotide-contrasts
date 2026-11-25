"""Layout-only drawing from the frozen, previously displayed source tables (MIT)."""
import ast
from pathlib import Path
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, Patch
import layout_common as lc

REPS={1:(lc.BLUE,'o'),2:(lc.GREEN,'s'),3:(lc.ORANGE,'^')}

def draw_f1(p):
    counts=pd.read_csv(p/'figure1_collection_scope.csv')
    hero=pd.read_csv(p/'figure1_qualified_group_counts.csv')
    fig=lc.newfig(height=4.05)
    a=lc.axes(fig,[.08,.61,.88,.29],letter='A');a.axis('off')
    labels=['Collected\nscore sets','Analyzable\nscore sets','Sets with\nmulti-key groups',
            'Missense\nscore rows','Historical\ngroup–key entries','Historical\nmulti-key groups']
    colors=[lc.BLUE,lc.GREEN,lc.ORANGE,lc.GREY,lc.BLUE,lc.GREEN]
    for i,row in counts.iterrows():
        x=(i%3)*.34;y=.87-(i//3)*.56
        a.text(x,y,f'{int(row.value):,}',fontsize=17,color=colors[i],ha='left',va='center',transform=a.transAxes)
        a.text(x,y-.17,labels[i],fontsize=10,ha='left',va='top',linespacing=1.12,transform=a.transAxes)
    b=lc.axes(fig,[.33,.15,.63,.30],letter='B',title='Groups meeting identity criteria')
    b.barh(range(3),hero.qualified_groups,color=colors[:3],height=.56)
    b.set_yticks(range(3),['Reference\ncompatible','Versioned\nrequest','Conditional\ntarget sequence']);b.invert_yaxis()
    b.set_xlim(0,6600);b.set_xticks([0,2000,4000,6000]);b.grid(axis='x',color=lc.LIGHT,lw=.6)
    b.set_xlabel('Qualified amino-acid groups',labelpad=6)
    for i,r in hero.iterrows():
        b.text(r.qualified_groups+100,i,f'{int(r.qualified_groups):,}  ({int(r.sets_with_qualified_groups)} sets)',va='center',fontsize=10,bbox={'facecolor':'white','edgecolor':'none','pad':1.2})
    fig.text(.33,.015,'Denominator: 13,041 historical multi-key groups',fontsize=9.5)
    lc.style(b)
    return fig

def draw_f2(p):
    g=pd.read_csv(p/'figure2_brca2_scatter_points.csv');assert len(g)==388 and g.n_variants.isin([2,3]).all()
    dr=pd.read_csv(p/'figure2_direction_interval_points.csv');ph=pd.read_csv(p/'figure2_phi_estimates.csv')
    fig=lc.newfig(height=5.3)
    a=lc.axes(fig,[.095,.23,.355,.66],letter='A',title='388 paired groups')
    truth=g.primary_concordant_bool.astype(str).str.lower().eq('true')
    for n,m in [(2,'o'),(3,'^')]:
        mask=g.n_variants.eq(2) if n==2 else g.n_variants.ge(3)
        a.scatter(g.loc[mask,'primary_d_Hap1_z'],g.loc[mask,'primary_d_mES_z'],s=23,
                  c=np.where(truth[mask],lc.BLUE,lc.ORANGE),marker=m,alpha=.8,lw=0)
    a.axhline(0,color=lc.LIGHT,lw=.9);a.axvline(0,color=lc.LIGHT,lw=.9)
    a.set_xlabel('HAP1 contrast\n(z-score units)');a.set_ylabel('mES contrast (z-score units)')
    hs=[Line2D([],[],marker=m,color='none',markerfacecolor=c,markeredgecolor=c,markersize=5,label=s)
        for m,c,s in [('o',lc.BLUE,'Concordant'),('o',lc.ORANGE,'Discordant'),('o',lc.GREY,'2 variants'),('^',lc.GREY,'3 variants')]]
    a.legend(handles=hs,loc='upper left',bbox_to_anchor=(-.05,-.20),ncol=2,frameon=False,
             handletextpad=.35,columnspacing=.65,labelspacing=.5,fontsize=9.5)
    b=lc.axes(fig,[.69,.59,.285,.30],letter='B',title='Concordance by pairing')
    methods=['primary_Hap1_extremes','reverse_mES_extremes','fixed_first_two_HGVS_keys']
    for j,method in enumerate(methods):
        for s,c,m,off in [('all',lc.BLUE,'o',-.20),('n2',lc.GREEN,'s',0),('n3plus',lc.ORANGE,'^',.20)]:
            r=dr[(dr.method==method)&(dr.stratum==s)].iloc[0];v=r.R_percent
            b.errorbar(v,j+off,xerr=[[v-r.ci_low_percent],[r.ci_high_percent-v]],fmt=m,color=c,ms=5,capsize=2,lw=1)
    b.set_yticks(range(3),['HAP1\nextremes','mES\nextremes','First two\nDNA keys']);b.set_ylim(2.5,-.5)
    b.set_xlim(30,70);b.set_xticks([30,50,70]);b.axvline(50,color=lc.GREY,ls='--',lw=.9)
    b.set_xlabel('Sign-concordant groups, R (%)',fontsize=9.5)
    hs=[Line2D([],[],marker=m,color=c,lw=0,markersize=5,label=s) for s,c,m in [('All',lc.BLUE,'o'),('2',lc.GREEN,'s'),('3',lc.ORANGE,'^')]]
    b.legend(handles=hs,loc='upper center',bbox_to_anchor=(.5,-.28),ncol=3,frameon=False,
             handletextpad=.25,columnspacing=.7,fontsize=9.5,title='DNA variants per group',title_fontsize=9.5)
    c=lc.axes(fig,[.69,.14,.285,.22],letter='C',title='Scale sensitivity')
    for i,r in ph.iterrows():
        col,m=REPS[i+1]
        if i==0:
            c.errorbar(r.phi_percent,i,xerr=[[r.phi_percent-r.ci_low_percent],[r.ci_high_percent-r.phi_percent]],fmt=m,color=col,ms=5,capsize=2,lw=1.1)
        else:c.scatter(r.phi_percent,i,s=30,color=col,marker=m)
        label=f'{r.phi_percent:.2f}' if i!=1 else f'{r.phi_percent:.3f}'
        c.text(r.phi_percent,i+.22,label,ha='center',va='top',fontsize=9.5)
    c.set_yticks(range(3),['Primary','Rank scale','Class control']);c.set_ylim(2.55,-.45)
    c.set_xlim(65,107);c.set_xticks([70,80,90,100]);c.set_xlabel('Phi (%)')
    fig.text(.59,.015,'Primary Phi 95% CI: 76.80–98.66%',fontsize=9.5)
    for ax in (a,b,c):lc.style(ax)
    return fig

def draw_f3(p):
    d=pd.read_csv(p/'figure3_msh2_all_endpoint_metrics.csv');assert len(d)==9
    fig=lc.newfig(height=3.75)
    a=lc.axes(fig,[.10,.24,.32,.61],letter='A',title='Primary endpoint')
    primary=d[d.endpoint.eq('6TG_mock')].sort_values('held_rep')
    for r in primary.itertuples():
        c,m=REPS[int(r.held_rep)];vs=[r.rmse_aa,r.rmse_full]
        a.plot([0,1],vs,color=c,marker=m,ms=5.5,lw=1.5)
        for x,v in enumerate(vs):
            dx=-7 if x==0 else 7
            # Outside endpoints, so adjacent replicate lines do not cover labels.
            a.annotate(f'{v:.2f}',(x,v),xytext=(dx,4 if int(r.held_rep)!=1 else 8),textcoords='offset points',
                       ha='right' if x==0 else 'left',va='bottom',fontsize=10,color=c)
    a.set_xlim(-.45,1.5);a.set_ylim(1.96,2.37);a.set_xticks([0,1],['Amino-acid\nmodel','DNA-informed\nmodel'])
    a.set_ylabel('6-TG / mock RMSE\n(log2 score units)');a.grid(axis='y',color=lc.LIGHT,lw=.6)
    b=lc.axes(fig,[.65,.24,.325,.61],letter='B',title='All endpoint gains')
    for j,e in enumerate(['6TG_mock','mock_P0','DB6_P0']):
        for r in d[d.endpoint.eq(e)].itertuples():
            c,m=REPS[int(r.held_rep)];v=r.relative_rmse_gain*100;y=j+{1:-.18,2:0,3:.18}[int(r.held_rep)]
            if e=='6TG_mock':
                lo,hi=[float(t)*100 for t in ast.literal_eval(r.relative_rmse_gain_95pct_ci)]
                b.errorbar(v,y,xerr=[[v-lo],[hi-v]],fmt=m,color=c,ms=5.5,capsize=2,lw=1.2)
            else:b.scatter(v,y,color=c,marker=m,s=30)
    b.axhspan(-.45,.45,color=lc.PALE_BLUE,zorder=0);b.axvline(0,color=lc.GREY,lw=.9,ls='--')
    b.set_yticks(range(3),['6-TG / mock','mock / P0','6-TG / P0']);b.set_ylim(2.5,-.5)
    b.set_xlim(-.3,7);b.set_xticks([0,2,4,6]);b.set_xlabel('Relative RMSE gain (%)')
    b.grid(axis='x',color=lc.LIGHT,lw=.6)
    fig.legend(handles=lc.legend_handles(),loc='lower center',bbox_to_anchor=(.55,.005),ncol=3,
               frameon=False,handletextpad=.4,columnspacing=1.4)
    for ax in (a,b):lc.style(ax)
    return fig

def draw_f4(p):
    null=pd.read_csv(p/'figure4_null_simulation_tau0.csv');real=pd.read_csv(p/'figure4_observed_workflows.csv')
    fig=lc.newfig(height=4.5)
    a=lc.axes(fig,[.10,.30,.35,.56],letter='A',title='Analytic two-variant rule')
    grid=np.linspace(-2.6,2.6,401);d1,d2=np.meshgrid(grid,grid)
    accepted=np.abs(d1)>1.96*np.sqrt(2*np.maximum((d1-d2)**2/4,1e-12))
    a.contourf(d1,d2,accepted.astype(int),levels=[.5,1.5],colors=[lc.PALE_GREEN])
    a.contour(d1,d2,accepted.astype(float),levels=[.5],colors=[lc.GREEN],linewidths=1)
    a.plot(grid,grid,color=lc.BLUE,lw=1,ls='--',label='D1 = D2')
    a.axhline(0,color=lc.GREY,lw=.7);a.axvline(0,color=lc.GREY,lw=.7)
    a.set_xlim(-2.6,2.6);a.set_ylim(-2.6,2.6);a.set_aspect('equal',adjustable='box')
    a.set_xlabel(r'Assay 1 contrast, $d_1$');a.set_ylabel(r'Assay 2 contrast, $d_2$')
    a.legend(handles=[Line2D([],[],color=lc.BLUE,lw=1,ls='--',label=r'$d_1 = d_2$'),
                      Patch(facecolor=lc.PALE_GREEN,edgecolor=lc.GREEN,label='Rule A acceptance')],
             loc='upper left',frameon=False,handlelength=1.2,fontsize=9.5)
    b=lc.axes(fig,[.72,.57,.255,.31],letter='B',title='Null simulation (τ/σ = 0)')
    cols=[lc.GREY,lc.BLUE,lc.GREEN,lc.ORANGE]
    for i,r in null.iterrows():
        b.plot([r.ci_low*100,r.ci_high*100],[i,i],color=cols[i],lw=2)
        b.scatter(r.mean_R*100,i,color=cols[i],s=32,zorder=3)
    b.set_yticks(range(4),[f'{s} (n̄={r.mean_selected:.1f})' for s,(_,r) in zip(['All','A','B','C'],null.iterrows())]);b.invert_yaxis()
    b.axvline(50,color=lc.GREY,ls='--',lw=.9);b.set_xlim(0,110);b.set_xticks([0,50,100])
    b.set_xlabel('Mean R (%)')
    c=lc.axes(fig,[.72,.14,.255,.27],letter='C',title='BRCA2 selection rules')
    labs=[]
    for i,r in real.iterrows():
        v=r.R*100;col=[lc.BLUE,lc.GREEN,lc.ORANGE,lc.GREY][i]
        if pd.notna(r.R_ci_low):c.errorbar(v,i,xerr=[[v-r.R_ci_low*100],[r.R_ci_high*100-v]],fmt='o',color=col,ms=5.5,capsize=2,lw=1.1)
        else:c.scatter(v,i,color=col,s=32)
        labs.append(['A','B-lab','B-tech','C'][i]+(f' (n̄={r.n_selected:.1f})' if i==3 else f' (n={int(r.n_selected)})'))
    c.set_yticks(range(4),labs);c.invert_yaxis();c.set_xlim(0,110);c.set_xticks([0,50,100])
    c.set_xlabel('Selected-group R (%)')
    for ax in (a,b,c):lc.style(ax)
    for ax in (b,c):ax.grid(axis='x',color=lc.LIGHT,lw=.6)
    return fig

if __name__=='__main__':
    inputs=Path(__file__).resolve().parent/'inputs'
    for fn,stem in [(draw_f1,'figure1_collection_testability'),(draw_f2,'figure2_brca2_paired_contrasts'),
                    (draw_f3,'figure3_msh2_controlled_endpoints'),(draw_f4,'figure4_structural_selection_counterexample')]:
        lc.finish(fn(inputs),stem);print('RENDERED',stem,flush=True)
