from babachi.bad_estimation import read_snps_file
from babachi.helpers import GenomeSNPsHandler
from babachi.visualize_segmentation import BabachiVisualizer
import seaborn as sns
from zipfile import ZipFile




snps, chrom_wrapper = read_snps_file(file_path='/home/ariuna/tcga_atacseq/processed_data/rssnps/BAM00225.snps.bed')

#self, snps_collection: GenomeSNPsHandler, BAD_file,

snps_collection = GenomeSNPsHandler(snps, chrom_wrapper)
visualizer = BabachiVisualizer(chromosomes_wrapper=chrom_wrapper)
visualizer.init_from_snps_collection(GenomeSNPsHandler,
                                             to_zip=False,
                                             ext='png', cosmic_file=None, cosmic_line=None,
                                             BAD_file='/home/ariuna/tcga_atacseq/babachi_all/BAM00225.snps.badmap.bed')