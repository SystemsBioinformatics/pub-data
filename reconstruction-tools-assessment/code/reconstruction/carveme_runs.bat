 

mkdir "D:\Dropbox\Research_Projects\Review_reconstruction\reconstructions\CA\lpl_run1"
cd D:\Dropbox\Research_Projects\Review_reconstruction\reconstructions\CA\lpl_run1
python D:\Programs\Anaconda2\Scripts\carve  "D:\Dropbox\Research_Projects\Review_reconstruction\reconstructions\Inputs\LPL\protein_fasta.faa" -o lp_001_ca.xml

mkdir "D:\Dropbox\Research_Projects\Review_reconstruction\reconstructions\CA\lpl_run2"
cd  D:\Dropbox\Research_Projects\Review_reconstruction\reconstructions\CA\lpl_run2
python D:\Programs\Anaconda2\Scripts\carve  "D:\Dropbox\Research_Projects\Review_reconstruction\reconstructions\Inputs\LPL\protein_fasta.faa" -u grampos -o lp_002_ca.xml



mkdir "D:\Dropbox\Research_Projects\Review_reconstruction\reconstructions\CA\bpe_run1"
cd D:\Dropbox\Research_Projects\Review_reconstruction\reconstructions\CA\bpe_run1
python D:\Programs\Anaconda2\Scripts\carve  "D:\Dropbox\Research_Projects\Review_reconstruction\reconstructions\Inputs\BPE\protein_fasta.faa" -o bp_001_ca.xml

mkdir "D:\Dropbox\Research_Projects\Review_reconstruction\reconstructions\CA\bpe_run2"
cd  D:\Dropbox\Research_Projects\Review_reconstruction\reconstructions\CA\bpe_run2
python D:\Programs\Anaconda2\Scripts\carve  "D:\Dropbox\Research_Projects\Review_reconstruction\reconstructions\Inputs\BPE\protein_fasta.faa" -u gramneg -o bp_002_ca.xml

cd D:\Dropbox\Research_Projects\Review_reconstruction\reconstructions\CA
