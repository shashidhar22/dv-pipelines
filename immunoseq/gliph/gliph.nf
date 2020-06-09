trb_files = Channel.fromPath(params.input.trb_path + '/*.tsv')
hla_files = Channel.fromPath(params.input.hla_path + "/*.tsv")
data_files = trb_files.merge(hla_files)

process gliph_analyze {
  echo false 

  publishDir "$params.output.folder/${tcr.getSimpleName()}", mode: 'copy'

  label 'gizmo_largenode'
    
  scratch "scratch"

  stageInMode 'copy'

  module "Perl"

  input:
    set path(tcr), path(hla) from data_files
    
  output:
    path "${tcr.getName()}-*.txt" into gliph_clonotypes
    
  """
  perl /home/sravisha/software/gliph/gliph/bin/gliph-group-discovery.pl --tcr ${tcr} 
  """

}