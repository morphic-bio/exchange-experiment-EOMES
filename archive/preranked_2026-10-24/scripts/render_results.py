"""One primary combined figure and one separately labeled ReMap sensitivity."""
import os
import argparse
import json
import tempfile
os.environ.setdefault('MPLCONFIGDIR',str(__import__('pathlib').Path(tempfile.gettempdir())/'eomes-release-mpl'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import ROOT, CONFIG, CONTRASTS, LABELS, TICKS, DE2_LABEL

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--single-degron',action='store_true',help='Additional primary figure with only the matched-control Degron comparison; retain original statistics.')
args=parser.parse_args()
suffix='_single_degron' if args.single_degron else ''
if args.single_degron:
    CONTRASTS=[c for c in CONTRASTS if c!='degron_parental']
    TICKS=['Whole-gene\nKO','Knock-in\nKO','Degron','CRISPRi']

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,
                     'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
sources=pd.read_csv(ROOT/'inputs/reference/signature_sources.tsv',sep='\t').set_index('signature')
sources.loc['TYSER_DE2','label']=DE2_LABEL
keys=CONFIG['cell_fate_display']

def target_heat(ax,data,label,show_y=True):
    values=data.pivot(index='direction',columns='contrast',values='fold_enrichment').reindex(index=['up','down'],columns=CONTRASTS)
    qs=data.pivot(index='direction',columns='contrast',values='BH_q').reindex(index=['up','down'],columns=CONTRASTS)
    im=ax.imshow(values,aspect='auto',cmap='YlGn',vmin=0,vmax=4)
    ax.set_yticks([0,1],['Upregulated DE','Downregulated DE'] if show_y else ['',''])
    ax.set_xticks(range(len(CONTRASTS)),TICKS,fontsize=9 if args.single_degron else 8)
    ax.tick_params(length=0)
    ax.set_title(label,loc='left',fontsize=11,fontweight='bold',pad=12)
    for i in range(2):
        for j in range(len(CONTRASTS)):
            val=values.iloc[i,j]
            ax.text(j,i,f'{val:.2f}'+('*' if qs.iloc[i,j]<.05 else ''),ha='center',va='center',fontsize=10,color='white' if val>2.8 else '#172030')
    return im

def fate_heat(ax,data,label,show_y=True):
    values=data.pivot(index='pathway',columns='contrast',values='NES').reindex(index=keys,columns=CONTRASTS)
    qs=data.pivot(index='pathway',columns='contrast',values='BH_q').reindex(index=keys,columns=CONTRASTS)
    im=ax.imshow(values,aspect='auto',cmap='RdBu_r',vmin=-5.5,vmax=5.5)
    ax.set_yticks(range(len(keys)),[sources.loc[k,'label'] for k in keys] if show_y else ['']*len(keys))
    ax.set_xticks(range(len(CONTRASTS)),TICKS,fontsize=9)
    ax.tick_params(length=0)
    ax.set_title(label,loc='left',fontsize=11,fontweight='bold',pad=12)
    for i in range(len(keys)):
        for j in range(len(CONTRASTS)):
            val=values.iloc[i,j]
            text='—' if pd.isna(val) else f'{val:.1f}'+('*' if qs.iloc[i,j]<.05 else '')
            ax.text(j,i,text,ha='center',va='center',fontsize=10 if args.single_degron else 8.5,color='white' if pd.notna(val) and abs(val)>3.3 else '#172030')
    return im

def save(fig,stem):
    for ext in ['png','pdf','svg']:
        fig.savefig(str(stem)+'.'+ext,dpi=260,bbox_inches='tight',pad_inches=.12,facecolor='white')
    plt.close(fig)

target=pd.read_csv(ROOT/'results/tables/chip48_directional_DE_enrichment.tsv',sep='\t')
fate=pd.read_csv(ROOT/'results/tables/cell_fate_GSEA.tsv',sep='\t')
fig,axes=plt.subplots(2,2,figsize=(12.5 if args.single_degron else 14.5,9.0),gridspec_kw={'height_ratios':[1,6.2],'hspace':.38,'wspace':.12})
fig.subplots_adjust(left=.29,right=.91,bottom=.105,top=.95)
im=target_heat(axes[0,0],target[target.top_n==500],'A   48 h ChIP targets: top 500')
target_heat(axes[0,1],target[target.top_n==1500],'48 h ChIP targets: top 1,500',False)
pos=axes[0,1].get_position()
cb=fig.colorbar(im,cax=fig.add_axes([.935,pos.y0,.014,pos.height]),ticks=[0,2,4]);cb.set_label('Fold enrichment',fontsize=9)
fate_heat(axes[1,0],fate[fate.variant=='all_genes'],'B   All expression-eligible genes')
im=fate_heat(axes[1,1],fate[fate.variant=='outside_ChIP48_top1500'],'C   Excluding 48 h ChIP top 1,500',False)
pos=axes[1,1].get_position()
cb=fig.colorbar(im,cax=fig.add_axes([.935,pos.y0+.04,.014,pos.height-.08]));cb.set_label('Normalized enrichment score (NES)')
save(fig,ROOT/f'results/figures/figure_1_chip48_and_cell_fate{suffix}')

if args.single_degron:
    caption=json.loads((ROOT/'results/figures/captions.json').read_text())['primary']
    caption['text']=caption['text'].replace(
        'both degron comparisons share perturbed samples, with matched uninduced or auxin-treated parental controls.',
        'Degron denotes the comparison with matched uninduced controls. The parental-control degron comparison is omitted from this display.')
    caption['text']+=' Enrichment values and adjusted p-values are retained from the original five-contrast analysis.'
    (ROOT/'results/figures/captions_single_degron.json').write_text(json.dumps({'primary':caption},indent=2)+'\n')
    (ROOT/'results/figures/figure_1_caption_single_degron.txt').write_text('Figure 1. '+caption['title']+' '+caption['text']+'\n')
    (ROOT/'qc/single_degron_figure.json').write_text(json.dumps(dict(
        displayed_contrasts=CONTRASTS,degron_control='matched uninduced',
        figure_size_inches=[12.5,9.0],upper_and_lower_value_font_pt=10,method_label_font_pt=9,
        omitted_contrast='degron_parental',statistics='Unchanged from original five-contrast analysis; original BH test families retained.'),indent=2)+'\n')
    print('Created additional primary figure with one Degron comparison.')
    raise SystemExit(0)

sens_target=pd.read_csv(ROOT/'sensitivity/remap/tables/directional_DE_enrichment.tsv',sep='\t')
sens_fate=pd.read_csv(ROOT/'sensitivity/remap/tables/cell_fate_GSEA.tsv',sep='\t')
fig,axes=plt.subplots(2,2,figsize=(14.5,9.0),gridspec_kw={'height_ratios':[1,6.2],'hspace':.38,'wspace':.12})
fig.subplots_adjust(left=.29,right=.91,bottom=.105,top=.95)
im=target_heat(axes[0,0],sens_target,'A   ReMap sensitivity: DE enrichment')
axes[0,1].axis('off')
pos=axes[0,0].get_position()
cb=fig.colorbar(im,cax=fig.add_axes([.935,pos.y0,.014,pos.height]),ticks=[0,2,4]);cb.set_label('Fold enrichment',fontsize=9)
fate_heat(axes[1,0],sens_fate[sens_fate.variant=='all_genes'],'B   All expression-eligible genes')
im=fate_heat(axes[1,1],sens_fate[sens_fate.variant=='outside_ReMap'],'C   Excluding ReMap targets',False)
pos=axes[1,1].get_position()
cb=fig.colorbar(im,cax=fig.add_axes([.935,pos.y0+.04,.014,pos.height-.08]));cb.set_label('Normalized enrichment score (NES)')
save(fig,ROOT/'sensitivity/remap/figures/figure_s1_remap_sensitivity')

caption_title='48-hour EOMES ChIP targets and cell-fate programs after EOMES perturbation.'
caption_text=(
    '(A) Enrichment of upregulated and downregulated DE genes among the fixed top 500 (left) and top 1,500 (right) '
    '48-hour ChIP targets. Values are observed/expected overlaps; * one-sided Fisher BH q < 0.05 across 20 tests. '
    '(B, C) Unweighted GSEA of cell-type signatures using all expression-eligible genes (B) or excluding the '
    '48-hour top 1,500 targets (C). Values are normalized enrichment scores (NES); red indicates higher and blue '
    'lower expression after perturbation. * BH q < 0.05 across all tested signatures and five contrasts, separately '
    'per panel. References: Tyser et al. (2021), human gastrula; Cao et al. (2020), human fetal atlas. DE2 denotes '
    'anterior-enriched definitive endoderm. Whole-gene and knock-in KO share WT samples; both degron comparisons '
    'share perturbed samples, with matched uninduced or auxin-treated parental controls. Reference similarity '
    'does not assign a cell fate; genes outside the target list are not proven indirect targets.'
)
captions=dict(primary=dict(title=caption_title,text=caption_text),sensitivity=dict(
    title='ReMap sensitivity analysis.',text='The primary workflow applied once to the fixed 1,446-gene ReMap EOMES reference set. Panel A uses directional Fisher tests, with BH correction across 10 comparisons; panel C excludes ReMap targets. Panel B is identical to the primary all-gene cell-fate analysis. This figure is supplementary and does not replace the 48-hour ChIP result. Other definitions and caveats follow Figure 1.'))
(ROOT/'results/figures/captions.json').write_text(json.dumps(captions,indent=2)+'\n')
(ROOT/'results/figures/figure_1_caption.txt').write_text('Figure 1. '+caption_title+' '+caption_text+'\n')
(ROOT/'sensitivity/remap/figures/figure_s1_caption.txt').write_text('Figure S1. '+captions['sensitivity']['title']+' '+captions['sensitivity']['text']+'\n')

# Gene-level support for the GSEA after excluding the 48h top-1500 targets.
genes=pd.read_csv(ROOT/'results/tables/gene_evidence.tsv.gz',sep='\t')
evidence=[]
for row in fate[fate.variant=='outside_ChIP48_top1500'].itertuples(index=False):
    if not isinstance(row.leadingEdge,str): continue
    for gene in row.leadingEdge.split(';'):
        evidence.append(dict(contrast=row.contrast,pathway=row.pathway,signature_label=sources.loc[row.pathway,'label'].replace('\n',' '),
                             signature_NES=row.NES,signature_BH_q=row.BH_q,symbol=gene))
ev=pd.DataFrame(evidence).merge(genes,on=['contrast','symbol'],how='left',validate='many_to_many')
assert not ev.ChIP48_top1500.any()
ev.to_csv(ROOT/'results/tables/outside_chip48_leading_edge_gene_evidence.tsv.gz',sep='\t',index=False,compression={'method':'gzip','mtime':0})
ev[(ev.signature_BH_q<.05)&ev.is_DE].to_csv(ROOT/'results/tables/outside_chip48_DE_program_genes.tsv',sep='\t',index=False)
print('Created the combined primary figure and one separate ReMap sensitivity figure.')
