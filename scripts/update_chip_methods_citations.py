"""Clarify peak calling and add software citations to the user's edited supplement."""
import argparse
import hashlib
import json
from pathlib import Path
from zipfile import ZipFile
from lxml import etree
from add_document_references import W, NS, sub, run, replace_in_paragraph

REPLACEMENTS={
    'Bowtie2 2.5.5 (default alignment settings)':'Bowtie2 2.5.5 [Langmead and Salzberg 2012] (default alignment settings)',
    'SAMtools 1.24 name-sorted':'SAMtools 1.24 [Danecek et al. 2021] name-sorted',
    'MACS2 2.2.9.1 called peaks':'MACS2 2.2.9.1 [Zhang et al. 2008] called peaks',
    'Thus, paired-end ChIP was treated as read-1 tags with fixed 200-bp extension, not BAMPE fragments.':
        'ChIP reads were aligned as paired-end data. For peak calling, MACS2 used the first mate of each pair and an assumed fragment length of 200 bp, rather than fragment lengths inferred from both mates.',
    'After blacklist exclusion, IDR 2.0.4.2 compared biological replicates':
        'After blacklist exclusion, the irreproducible discovery rate (IDR) method [Li et al. 2011], implemented in IDR 2.0.4.2, compared biological replicates',
}
REFERENCES=[
    ('[Danecek et al. 2021] ', 'Danecek, P., Bonfield, J. K., Liddle, J., et al. (2021). Twelve years of SAMtools and BCFtools. GigaScience, 10(2), giab008. ', '10.1093/gigascience/giab008'),
    ('[Langmead and Salzberg 2012] ', 'Langmead, B., and Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie 2. Nature Methods, 9, 357–359. ', '10.1038/nmeth.1923'),
    ('[Li et al. 2011] ', 'Li, Q., Brown, J. B., Huang, H., and Bickel, P. J. (2011). Measuring reproducibility of high-throughput experiments. The Annals of Applied Statistics, 5(3), 1752–1779. ', '10.1214/11-AOAS466'),
    ('[Zhang et al. 2008] ', 'Zhang, Y., Liu, T., Meyer, C. A., et al. (2008). Model-based Analysis of ChIP-Seq (MACS). Genome Biology, 9, R137. ', '10.1186/gb-2008-9-9-r137'),
]

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('input',type=Path);ap.add_argument('output',type=Path)
    args=ap.parse_args();assert args.input.resolve()!=args.output.resolve()
    with ZipFile(args.input) as z:
        root=etree.fromstring(z.read('word/document.xml'));body=root.find('w:body',NS)
        children=list(body);before=[etree.tostring(e) for e in children]
        changes=[]
        for i,p in enumerate(children):
            if p.tag!='{'+W+'}p':continue
            for old,new in REPLACEMENTS.items():
                if replace_in_paragraph(p,old,new):changes.append(dict(paragraph=i,before=old,after=new))
        assert len(changes)==5,changes
        for i,p in enumerate(children):
            if i not in {c['paragraph'] for c in changes}:assert etree.tostring(p)==before[i]
        text=lambda p:''.join(p.xpath('.//w:t/text()',namespaces=NS))
        heading=next(p for p in children if text(p)=='References')
        existing=[p for p in children[children.index(heading)+1:] if p.tag=='{'+W+'}p' and text(p).startswith('[')]
        assert len(existing)==4
        relroot=etree.fromstring(z.read('word/_rels/document.xml.rels'))
        relns='http://schemas.openxmlformats.org/package/2006/relationships'
        next_id=max(int(e.get('Id')[3:]) for e in relroot if e.get('Id','').startswith('rId'))+1
        for label,entry,doi in REFERENCES:
            p=etree.Element('{'+W+'}p');pr=sub(p,'pPr')
            sub(pr,'pStyle',val='Methodsbody');sub(pr,'keepLines');sub(pr,'spacing',after='180')
            run(p,label,bold=True);run(p,entry);sub(sub(p,'r'),'br')
            relid=f'rId{next_id}';next_id+=1
            etree.SubElement(relroot,'{'+relns+'}Relationship',Id=relid,
                Type='http://schemas.openxmlformats.org/officeDocument/2006/relationships/hyperlink',
                Target='https://doi.org/'+doi,TargetMode='External')
            h=sub(p,'hyperlink');h.set('{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id',relid)
            run(h,'https://doi.org/'+doi)
            successors=[x for x in existing if text(x).split(']')[0]>label.strip().split(']')[0]]
            before_node=successors[0] if successors else body.find('w:sectPr',NS)
            body.insert(body.index(before_node),p)
            existing=sorted(existing+[p],key=text)
        modified={'word/document.xml':etree.tostring(root,xml_declaration=True,encoding='UTF-8',standalone=True),
                  'word/_rels/document.xml.rels':etree.tostring(relroot,xml_declaration=True,encoding='UTF-8',standalone=True)}
        with ZipFile(args.output,'w') as out:
            for info in z.infolist():out.writestr(info,modified.get(info.filename,z.read(info.filename)))
    with ZipFile(args.input) as a,ZipFile(args.output) as b:
        assert all(a.read(n)==b.read(n) for n in a.namelist() if n not in modified)
    audit=dict(source=str(args.input),output=str(args.output),source_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),
        output_sha256=hashlib.sha256(args.output.read_bytes()).hexdigest(),changes=changes,reference_count=8,
        all_other_paragraphs_and_existing_references_preserved=True,images_styles_and_section_layout_preserved=True)
    (Path(__file__).resolve().parents[1]/'qc/chip_methods_citation_edit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(args.output)

if __name__=='__main__':main()
