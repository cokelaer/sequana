from sequana import bedtools
from sequana.utils.datatables_js import DataTable, DataTableFunction

from . import test_dir

def test_datatables():
        bed = bedtools.GenomeCov(f"{test_dir}/data/bed/JB409847.bed",
                                 f"{test_dir}/data/genbank/JB409847.gbk")
        fasta = f"{test_dir}/data/fasta/JB409847.fasta"
        bed.compute_gc_content(fasta)

        c = bed.chr_list[0]
        c.run(4001)
        rois = c.get_rois()
        rois.df['link'] = 'test'
        datatable_js = DataTableFunction(rois.df, 'roi')
        datatable_js.set_links_to_column('link', 'start')
        datatable_js.datatable_options = {'scrollX': 'true',
                                          'pageLength': 15,
                                          'scrollCollapse' : 'true',
                                          'dom': 'Bfrtip',
                                          'buttons': ['copy', 'csv']}
        datatable = DataTable(rois.df, 'rois', datatable_js)
        html_table = datatable.create_datatable(float_format='%.3g')
