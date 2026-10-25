"""Add references to an edited Word document without rebuilding its content."""
import argparse
import hashlib
import json
from pathlib import Path
from zipfile import ZipFile
from lxml import etree

W='http://schemas.openxmlformats.org/wordprocessingml/2006/main'
XML='http://www.w3.org/XML/1998/namespace'
NS={'w':W}

REPLACEMENTS={
    'Luo et al., Nature Genetics (2023)':'[Luo et al. 2023]',
    'Li et al., Nature Genetics (2019)':'[Li et al. 2019]',
    'Tyser et al., Nature (2021)':'[Tyser et al. 2021]',
    'Cao et al., Science (2020)':'[Cao et al. 2020]',
    'Tyser et al. (2021)':'[Tyser et al. 2021]',
    'Cao et al. (2020)':'[Cao et al. 2020]',
}
REFERENCES=[
    ('[Cao et al. 2020] ', "Cao, J., O’Day, D. R., Pliner, H. A., et al. (2020). A human cell atlas of fetal gene expression. Science, 370(6518), eaba7721. ", '10.1126/science.aba7721'),
    ('[Li et al. 2019] ', 'Li, Q. V., Dixon, G., Verma, N., et al. (2019). Genome-scale screens identify JNK–JUN signaling as a barrier for pluripotency exit and endoderm differentiation. Nature Genetics, 51, 999–1010. ', '10.1038/s41588-019-0408-9'),
    ('[Luo et al. 2023] ', 'Luo, R., Yan, J., Oh, J. W., et al. (2023). Dynamic network-guided CRISPRi screen identifies CTCF-loop-constrained nonlinear enhancer gene regulatory activity during cell state transitions. Nature Genetics, 55, 1336–1346. ', '10.1038/s41588-023-01450-7'),
    ('[Tyser et al. 2021] ', 'Tyser, R. C. V., Mahammadov, E., Nakanoh, S., et al. (2021). Single-cell transcriptomic characterization of a gastrulating human embryo. Nature, 600, 285–289. ', '10.1038/s41586-021-04158-y'),
]

def sub(parent,tag,**attrs):
    e=etree.SubElement(parent,'{'+W+'}'+tag)
    for k,v in attrs.items():e.set('{'+W+'}'+k,str(v))
    return e

def run(parent,text,bold=False):
    r=sub(parent,'r')
    if bold:sub(sub(r,'rPr'),'b')
    t=sub(r,'t');t.set('{'+XML+'}space','preserve');t.text=text

def replace_in_paragraph(p,old,new):
    nodes=p.findall('.//w:t',NS)
    full=''.join(n.text or '' for n in nodes)
    start=full.find(old)
    if start<0:return False
    assert full.count(old)==1
    end=start+len(old);pos=0;inserted=False
    for n in nodes:
        s=n.text or '';nend=pos+len(s)
        if pos<end and nend>start:
            left=max(0,start-pos);right=min(len(s),end-pos)
            n.text=s[:left]+(new if not inserted else '')+s[right:]
            n.set('{'+XML+'}space','preserve');inserted=True
        pos=nend
    assert ''.join(n.text or '' for n in nodes)==full.replace(old,new)
    return True

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('input',type=Path);ap.add_argument('output',type=Path)
    args=ap.parse_args()
    assert args.input.resolve()!=args.output.resolve(), 'Preserve the user-edited source.'
    with ZipFile(args.input) as z:
        root=etree.fromstring(z.read('word/document.xml'))
        body=root.find('w:body',NS)
        original_children=list(body)
        before=[etree.tostring(e) for e in original_children]
        changes=[]
        for i,p in enumerate(original_children):
            if p.tag!='{'+W+'}p':continue
            for old,new in REPLACEMENTS.items():
                if replace_in_paragraph(p,old,new):changes.append({'paragraph':i,'before':old,'after':new})
        assert len(changes)==6, changes
        # All non-citation paragraphs and section settings remain byte-identical in XML.
        changed_indices={x['paragraph'] for x in changes}
        for i,e in enumerate(original_children):
            if i not in changed_indices:assert etree.tostring(e)==before[i]
        sect=body.find('w:sectPr',NS)
        insertion=body.index(sect)
        heading=etree.Element('{'+W+'}p')
        props=sub(heading,'pPr');sub(props,'pStyle',val='Methodstitle')
        sub(props,'keepNext');sub(props,'pageBreakBefore')
        run(heading,'References')
        body.insert(insertion,heading);insertion+=1
        for label,text,doi in REFERENCES:
            p=etree.Element('{'+W+'}p')
            props=sub(p,'pPr');sub(props,'pStyle',val='Methodsbody')
            sub(props,'keepLines');sub(props,'spacing',after='180')
            run(p,label,bold=True);run(p,text)
            sub(sub(p,'r'),'br')
            # Existing DOI relationships are reused; no package parts are regenerated.
            relroot=etree.fromstring(z.read('word/_rels/document.xml.rels'))
            relid=next(e.get('Id') for e in relroot if e.get('Target')=='https://doi.org/'+doi)
            h=sub(p,'hyperlink');h.set('{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id',relid)
            run(h,'https://doi.org/'+doi)
            body.insert(insertion,p);insertion+=1
        newxml=etree.tostring(root,xml_declaration=True,encoding='UTF-8',standalone=True)
        with ZipFile(args.output,'w') as out:
            for info in z.infolist():out.writestr(info,newxml if info.filename=='word/document.xml' else z.read(info.filename))
    with ZipFile(args.input) as a,ZipFile(args.output) as b:
        assert a.namelist()==b.namelist()
        assert all(a.read(n)==b.read(n) for n in a.namelist() if n!='word/document.xml')
    audit={'source':str(args.input),'source_sha256':hashlib.sha256(args.input.read_bytes()).hexdigest(),
           'output':str(args.output),'output_sha256':hashlib.sha256(args.output.read_bytes()).hexdigest(),
           'citation_replacements':changes,'reference_count':4,'noncitation_content_and_formatting_preserved':True,
           'all_non_document_xml_package_parts_identical':True}
    qc=Path(__file__).resolve().parents[1]/'qc/document_references_edit.json'
    qc.write_text(json.dumps(audit,indent=2)+'\n')
    print(args.output)

if __name__=='__main__':main()
