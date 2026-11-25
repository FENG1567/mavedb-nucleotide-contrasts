"""Shared print typography and layout controls (original code, MIT)."""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.text import Text

ROOT = Path(__file__).resolve().parent
OUTPUT_DIR = ROOT / 'figures'
for font in (ROOT / 'fonts').glob('*.ttf'):
    font_manager.fontManager.addfont(str(font))
BLUE='#3B6F9C'; GREEN='#5B8E75'; ORANGE='#C7864D'; GREY='#414950'; LIGHT='#D9E0E5'
INK='#17232D'; PALE_BLUE='#E7EEF4'; PALE_GREEN='#E7EFEA'
mpl.rcParams.update({'font.family':'Arial','font.weight':'bold','font.size':10,
 'axes.labelsize':10,'axes.labelweight':'bold','axes.titleweight':'bold','axes.titlesize':11,
 'xtick.labelsize':9.5,'ytick.labelsize':9.5,'legend.fontsize':9.5,
 'text.color':INK,'axes.labelcolor':INK,'xtick.color':INK,'ytick.color':INK,
 'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.85,
 'xtick.major.width':.85,'ytick.major.width':.85,'svg.fonttype':'none','pdf.fonttype':42,
 'mathtext.fontset':'custom','mathtext.rm':'Arial:bold','mathtext.it':'Arial:bold',
 'mathtext.bf':'Arial:bold','figure.facecolor':'white','savefig.facecolor':'white'})

def newfig(width=7.2,height=4):
    return plt.figure(figsize=(width,height))

def axes(fig,rect,letter=None,title=None):
    ax=fig.add_axes(rect)
    if letter: tag(fig,ax,letter)
    if title: ax.set_title(title,loc='left',pad=10)
    return ax

def tag(fig,ax,letter):
    b=ax.get_position()
    # Fixed print-space type, independent of the size of the individual panel.
    t=fig.text(max(.012,b.x0-.068),b.y1+.018,letter.upper(),fontsize=18,
               fontname='Arial',fontweight='bold',ha='left',va='bottom')
    t.set_gid('panel-label-'+letter.upper())
    return t

def style(ax):
    ax.tick_params(length=3,pad=3)
    ax.spines['left'].set_color(GREY);ax.spines['bottom'].set_color(GREY)
    ax.set_axisbelow(True)
    for t in ax.get_xticklabels()+ax.get_yticklabels(): t.set_fontweight('bold')

def legend_handles():
    from matplotlib.lines import Line2D
    return [Line2D([],[],color=c,marker=m,lw=0,markersize=5,label=f'Held-out {i}')
            for i,c,m in zip((1,2,3),(BLUE,GREEN,ORANGE),('o','s','^'))]

def finish(fig,stem):
    out=OUTPUT_DIR;out.mkdir(parents=True,exist_ok=True)
    texts=[]
    for t in fig.findobj(Text):
        if not t.get_text(): continue
        t.set_fontweight('bold');t.set_fontfamily('Arial')
        texts.append({'text':t.get_text(),'size_pt':t.get_fontsize(),'weight':t.get_fontweight(),'gid':t.get_gid()})
    fig.canvas.draw()
    # Do not use tight bounding-box export: preserve exact panel-label point size.
    fig.savefig(out/(stem+'.pdf'),metadata={'Creator':'Matplotlib','Title':stem})
    fig.savefig(out/(stem+'.svg'))
    fig.savefig(out/(stem+'.png'),dpi=300)
    fig.savefig(out/(stem+'_embed.png'),dpi=450)
    (out/(stem+'_typography.json')).write_text(json.dumps(texts,indent=2),encoding='utf8')
    plt.close(fig)
