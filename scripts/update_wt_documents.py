"""Update only GSEA methods, captions and figure images in the archived documents."""
import copy
import json
from pathlib import Path
from zipfile import ZipFile
from lxml import etree
from common import ROOT, sha256

W = 'http://schemas.openxmlformats.org/wordprocessingml/2006/main'
NS = {'w': W, 'a': 'http://schemas.openxmlformats.org/drawingml/2006/main',
      'wp': 'http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing',
      'r': 'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}

METHODS = [
    'For each contrast, all expression-eligible genes were ranked by decreasing signed log FC. Eligibility required q75_pass and max(mean KO CPM, mean control CPM) ≥ 1; DE significance was not required. For repeated symbols, the eligible locus with highest AveExpr was retained, with Ensembl ID breaking ties. Unweighted GSEA (gseaParam = 0) tested signatures with 10–1,000 represented genes. These observed rankings were preserved and WT coefficients were matched by Ensembl gene ID.',
    'The empirical reference used the nine selected EOMES WT libraries and 1,566 usable WT-versus-WT 3v3 assignments from GSEA_expanded_WT_3v3. Each six-library subset was TMM-normalized on the original 23,100-gene q75 universe; limma-voom fitted source group plus pseudo-condition. There were 783 fitted assignments, scored in both label orientations after reranking. Each contrast’s observed expression mask was fixed across all assignments. Reference reconstruction reproduced all 7,830 saved upstream custom-ChIP scores.',
    'For observed enrichment score ES, the p-value was (1 + number of WT scores at least as extreme)/(1 + number of same-sign WT scores). NES was ES divided by the mean absolute same-sign WT score. Positive/negative NES indicates enrichment toward increased/decreased expression. For target-excluded GSEA, the top-1,500 ChIP symbols were removed from both observed and WT rankings before scoring. BH correction included all tested signatures and all five contrasts, separately for all genes (441 tests) and target exclusion (436 tests); figure stars use BH_q < 0.05.',
    'Random-gene-set fgseaMultilevel results and balanced-only/residual-df ≥ 2 WT sensitivities are retained in the tables. The WT reference is empirical, not permutation of actual KO/control labels: WT samples were PCA-selected, assignments overlap, and six-library reference variance may differ from full observed comparisons. Enrichment of a marker program does not establish a cell fate.'
]


def set_text(p, text):
    nodes = p.findall('.//w:t', NS)
    assert nodes
    nodes[0].text = text
    nodes[0].set('{http://www.w3.org/XML/1998/namespace}space','preserve')
    for node in nodes[1:]:
        node.text = ''


def main():
    audits=[]
    for name, suffix in [('EOMES_supp_updated.docx','_single_degron'),('figure_1_with_caption.docx','')]:
        source=ROOT/'archive/preranked_2026-10-24/results/documents'/name
        target=ROOT/'results/documents'/name
        caption=json.loads((ROOT/f'results/figures/captions{suffix}.json').read_text())['primary']
        caption_text='Figure 1. '+caption['title']+' '+caption['text']
        if suffix:
            caption_text=caption_text.replace('Tyser et al. (2021)','[Tyser et al. 2021]').replace('Cao et al. (2020)','[Cao et al. 2020]')
        with ZipFile(source) as zin:
            doc=etree.fromstring(zin.read('word/document.xml'))
            body=doc.find('w:body',NS)
            paragraphs=body.findall('w:p',NS)
            changed=[]
            untouched=[]
            for i,p in enumerate(paragraphs):
                text=''.join(t.text or '' for t in p.findall('.//w:t',NS))
                if not suffix and p.findall('.//a:blip',NS):
                    extent=p.find('.//wp:extent',NS)
                    ratio=(8.2*914400)/int(extent.get('cx'))
                    for e in p.findall('.//wp:extent',NS)+p.findall('.//a:xfrm/a:ext',NS):
                        for dim in ['cx','cy']:e.set(dim,str(round(int(e.get(dim))*ratio)))
                    changed.append(i)
                elif text.startswith('Figure 1.'):
                    # Preserve bold title and regular caption body.
                    nodes=p.findall('.//w:t',NS)
                    assert len(nodes)>=2
                    nodes[0].text='Figure 1. '+caption['title']+' '
                    nodes[1].text=caption_text[len(nodes[0].text):]
                    for t in nodes[2:]:t.text=''
                    changed.append(i)
                elif text.startswith('For each contrast, all expression-eligible genes were ranked'):
                    set_text(p,METHODS[0]);at=body.index(p)+1
                    for text in METHODS[1:]:
                        extra=copy.deepcopy(p);set_text(extra,text);body.insert(at,extra);at+=1
                    changed.append(i)
                elif text in ['Cell-type program enrichment by DE signals','Reference']:
                    props=p.find('w:pPr',NS)
                    if props is None:
                        props=etree.Element('{'+W+'}pPr');p.insert(0,props)
                    if props.find('w:keepNext',NS) is None:etree.SubElement(props,'{'+W+'}keepNext')
                    if text=='Reference' and props.find('w:pageBreakBefore',NS) is None:
                        etree.SubElement(props,'{'+W+'}pageBreakBefore')
                    changed.append(i)
                else:
                    untouched.append((p,etree.tostring(p)))
            blips=doc.findall('.//a:blip',NS);assert len(blips)==1
            relid=blips[0].get('{'+NS['r']+'}embed')
            rels=etree.fromstring(zin.read('word/_rels/document.xml.rels'))
            media='word/'+next(r.get('Target') for r in rels if r.get('Id')==relid)
            assert all(etree.tostring(p)==before for p,before in untouched)
            replacements={'word/document.xml':etree.tostring(doc,xml_declaration=True,encoding='UTF-8',standalone=True),
                          media:(ROOT/f'results/figures/figure_1_chip48_and_cell_fate{suffix}.png').read_bytes()}
            with ZipFile(target,'w') as zout:
                for item in zin.infolist():zout.writestr(item,replacements.get(item.filename,zin.read(item.filename)))
        audits.append(dict(source=str(source.relative_to(ROOT)),source_sha256=sha256(source),
                           output=str(target.relative_to(ROOT)),output_sha256=sha256(target),changed_original_paragraphs=changed,
                           references_and_other_paragraphs_preserved=True))
    (ROOT/'qc/wt_document_update.json').write_text(json.dumps(audits,indent=2)+'\n')


if __name__=='__main__':main()
