#### Congfiguration ####

# fetch the data from the paper, parse it into something sensibke

require 'catpaws/ec2'
require 'pathname'

set :aws_access_key,  ENV['AMAZON_ACCESS_KEY']
set :aws_secret_access_key , ENV['AMAZON_SECRET_ACCESS_KEY']
set :ec2_url, ENV['EC2_URL']
set :key, ENV['EC2_KEY'] #this should be the name of the key in aws, not the actual keyfile
set :ssh_options, {
  :keys => [ENV['EC2_KEYFILE']],
  :user => "ubuntu"
}
set :ami, 'ami-52794c26' #32-bit ubuntu lucid server (eu-west-1)
set :instance_type, 'm1.small'
set :working_dir, '/mnt/work'

set :group_name, 'Johnson_ChIPPET'
set :nhosts, 1

set :snap_id, `cat SNAPID`.chomp #empty until you've created a snapshot
set :vol_id, `cat VOLUMEID`.chomp #empty until you've created a new volume
set :ebs_size, 3  
set :availability_zone, 'eu-west-1a'  #wherever your ami is. 
set :dev, '/dev/sdf'
set :mount_point, '/mnt/data'

set :git_url, 'http://github.com/cassj/Johnson_CHIPPET/raw/master'

# Try and load a local config file to override any of the above values, should one exist.
# So that if you change these values, they don't get overwritten if you update the repos.
begin
 load("Capfile.local")
rescue Exception
end

#cap EC2:start
#cap EBS:create (unless you want to use the one that already exists)
#cap EBS:attach
#cap EBS:format_xfs (unless you're using a snapshot)
#cap EBS:mount_xfs
#
# run tasks as required
# 
#cap EBS:snapshot
#cap EBS:unmount
#cap EBS:detach
#cap EBS:delete
#cap EC2:stop


#### Tasks ####


desc "Fetch data from publication"
task :publication_data, :roles => group_name do
  
  run "mkdir -p #{mount_point}/publication"
  run "cd #{mount_point}/publication && curl 'http://www.plosbiology.org/article/fetchObjectAttachment.action?uri=info%3Adoi%2F10.1371%2Fjournal.pbio.0060256&representation=PDF' > Johnson.pdf"
  run "cd #{mount_point}/publication && curl 'http://www.plosbiology.org/article/fetchSingleRepresentation.action?uri=info:doi/10.1371/journal.pbio.0060256.sd007' >  SupplementalMethods.txt"
  for i in 1..9
    run "cd #{mount_point}/publication && curl 'http://www.plosbiology.org/article/fetchSingleRepresentation.action?uri=info:doi/10.1371/journal.pbio.0060256.sg00#{i}' > FigS#{i} "
  end

  for i in 10..13
    run "cd #{mount_point}/publication && curl 'http://www.plosbiology.org/article/fetchSingleRepresentation.action?uri=info:doi/10.1371/journal.pbio.0060256.sg0#{i}' >FigS#{i} "
  end
  for i in 1..6
    run "cd #{mount_point}/publication && curl 'http://www.plosbiology.org/article/fetchSingleRepresentation.action?uri=info:doi/10.1371/journal.pbio.0060256.sd00#{i}' > DatasetS#{i}"
  end  
  put("Figure S1. Validation of ESC Pluripotency by Cell Sorting\nFigure S2. Expression of Neural Stem Cell Markers by NS5 Cells\nFigure S3. PCR Validation of ChIP-chip in ESC\nFigure S4. qPCR Validation of ChIP-chip Data: Shared ESC/NSC RE1s\nFigure S5. qPCR Validation of ChIP-chip Data: ESC-Specific Binding Sites\nFigure S6. qPCR Validation of ChIP-chip Data: NSC-Specific RE1s\nFigure S7. qPCR Validation of ChIP-chip Data: 3T3-Specific RE1s\nFigure S8. qPCR Validation of ChIP-PET to Identify PET Cluster Size Cutoff (REST ChIP in ESC)\nFigure S9. Comparison of ChIP-PET and ChIP-chip REST Binding Predictions\nFigure S10. qPCR Validation of ChIP-chip Data: 'No Motif' Binding Sites in ESC\nFigure S11. qPCR Validation of ChIP-PET Data: Shared ESC/NSC RE1s\nFigure S12. qPCR Validation of ChIP-PET Data: ESC-Specific Binding Sites\nFigure S13. REST Represses Most Effectively from Promoter-Proximal Binding Sites\nFigure S14. Regulation of Pluripotency Genes by REST\nFigure S15. No Evidence for REST Recruitment to the Mouse mir-21 Locus\n\nDataset S1. ChIP-chip Data\nDataset S2. ESC ChIP-PET Motifs\nDataset S3. NSC ChIP-PET Motifs\nDataset S4. ChIP-PET Target Genes\nDataset S5. DN:REST Gene Expression Data\nDataset S6. ESC-Specific Target Genes", "#{mount_point}/publication/readme")
end
before 'publication_data','EC2:start'


desc "install R on all running instances in group group_name"
task :install_r, :roles  => group_name do
  user = variables[:ssh_options][:user]
  sudo 'apt-get update'
  sudo 'apt-get -y install r-base'
  sudo 'apt-get -y install build-essential libxml2 libxml2-dev libcurl3 libcurl4-openssl-dev'
  run "wget --no-check-certificate '#{git_url}/scripts/R_setup.R' -O #{working_dir}/R_setup.R"
  run "cd #{working_dir} && chmod +x R_setup.R"
  sudo "Rscript #{working_dir}/R_setup.R"
end
before "install_r", "EC2:start"


### Process the expression data ###

set :esc_xpn_file, "ESC-REST-DN-raw.csv"
set :ns5_xpn_file, "NS5-REST-DN-raw.csv"


desc "Fetch expression data from Buckley Lab server"
task :expression_data, :roles => group_name do  
  run "mkdir -p #{mount_point}/publication"
  upload("data/E14-REST_DN_QC_List-w-Replicates.csv", "#{mount_point}/publication/#{esc_xpn_file}")
  upload("data/NS5-REST-DN-raw.csv", "#{mount_point}/publication/#{ns5_xpn_file}")
  #create results dirs
  run "mkdir -p #{mount_point}/NS5/PET"
  run "mkdir -p #{mount_point}/ESC/PET"
  run "mkdir -p #{mount_point}/NS5/expression"
  run "mkdir -p #{mount_point}/ESC/expression"

end
before 'expression_data', 'EC2:start'



desc "run QC checks on the raw data"
task :qc_expression_data, :roles => group_name do
  puts "TODO"
end
before "qc_expression_data", "EC2:start"


desc "run pre-processing on expression data"
task :pp_expression_data, :roles => group_name do
  run "mkdir -p #{working_dir}/scripts"
  run "cd #{working_dir}/scripts && curl #{git_url}/scripts/limma_xpn.R > limma_xpn.R"
  run "chmod +x #{working_dir}/scripts/limma_xpn.R"
  run "cd #{mount_point}/NS5/expression && Rscript #{working_dir}/scripts/limma_xpn.R #{mount_point}/publication/#{ns5_xpn_file} limma_results.csv"
  run "cd #{mount_point}/ESC/expression && Rscript #{working_dir}/scripts/limma_xpn.R #{mount_point}/publication/#{esc_xpn_file} limma_results.csv"
end
before "pp_expression_data", "EC2:start"


desc "run QC checks on the pre-processed quality control"
task "pp_qc_expression_data", :roles => group_name do
  puts "TODO"
end  
before "pp_qc_exprssion_data", "EC2:start"


desc "Fetch ReMoat data which has mm9 probe positions"
task :get_remoat_anno, :roles => group_name do
  run "mkdir -p #{working_dir}/lib"
  run "rm -Rf  #{working_dir}/lib/Annotation_Illumina_Mouse*"
  run "cd #{working_dir}/lib && curl http://www.compbio.group.cam.ac.uk/Resources/Annotation/final/Annotation_Illumina_Mouse-WG-V1_mm9_V1.0.0_Aug09.zip > Annotation_Illumina_Mouse-WG-V1_mm9_V1.0.0_Aug09.zip "
  run "cd #{working_dir}/lib && unzip Annotation_Illumina_Mouse-WG-V1_mm9_V1.0.0_Aug09.zip"
end 
before 'get_remoat_anno', 'EC2:start'


desc "Make an IRanges RangedData object from expression data"
task :xpn2rd, :roles => group_name do
  user = variables[:ssh_options][:user]
  run "cd #{working_dir}/scripts && curl '#{git_url}/scripts/xpn_csv_to_iranges.R' >  xpn_csv_to_iranges.R"
  run "cd #{working_dir}/scripts && chmod +x xpn_csv_to_iranges.R"
  run "cd #{mount_point}/ESC/expression && Rscript #{working_dir}/scripts/xpn_csv_to_iranges.R limma_results.csv #{working_dir}/lib/Annotation_Illumina_Mouse-WG-V1_mm9_V1.0.0_Aug09.txt"
  run "cd #{mount_point}/NS5/expression && Rscript #{working_dir}/scripts/xpn_csv_to_iranges.R limma_results.csv #{working_dir}/lib/Annotation_Illumina_Mouse-WG-V1_mm9_V1.0.0_Aug09.txt"

end
before "xpn2rd","EC2:start"  


desc "Fetch expression data results"
task :get_xpn, :roles=> group_name do
  `mkdir -p results/NS5/expression`
  `mkdir -p results/ESC/expression`
  download("#{mount_point}/ESC/expression/limma_rd.csv", "results/ESC/expression/limma_rd.csv")
  download("#{mount_point}/NS5/expression/limma_rd.csv", "results/NS5/expression/limma_rd.csv")

end 
before "get_xpn", "EC2:start"



### Process the PET data ###

desc "install liftOver"
task :install_liftover, :roles => group_name do
  run "curl http://hgdownload.cse.ucsc.edu/admin/exe/linux.i386/liftOver > #{working_dir}/liftOver"
  run "sudo mv #{working_dir}/liftOver /usr/bin/liftOver"
  run "sudo chmod +x /usr/bin/liftOver"
  run "mkdir -p #{working_dir}/lib"
  run "curl http://hgdownload.cse.ucsc.edu/goldenPath/mm8/liftOver/mm8ToMm9.over.chain.gz > #{working_dir}/lib/mm8ToMm9.over.chain.gz"
  run "cd #{working_dir}/lib && gunzip mm8ToMm9.over.chain.gz"
  
end
before "install_liftover", "EC2:start"



desc "Parse the PET data to csv"
task :parse_pet, :roles => group_name do
  run "mkdir -p #{working_dir}/scripts"
  run "mkdir -p #{mount_point}/ESC/PET"
  run "mkdir -p #{mount_point}/NS5/PET"
  run "cd #{working_dir}/scripts && curl #{git_url}/scripts/convert_excel.pl > convert_excel.pl"
  run "cd #{working_dir}/scripts && chmod +x convert_excel.pl"
  run "sudo apt-get -y install  libspreadsheet-parseexcel-perl"
  run "cd #{mount_point}/publication &&  #{working_dir}/scripts/convert_excel.pl DatasetS2 #{mount_point}/ESC/PET"
  run "cd #{mount_point}/publication &&  #{working_dir}/scripts/convert_excel.pl DatasetS3 #{mount_point}/NS5/PET"
end
before 'parse_pet', 'EC2:start'


# Merge the PET csv files from different RE1 site types and convert to IRanges
desc "PET files to mm9 Iranges"
task :pet2iranges, :roles => group_name do
  run "mkdir -p #{working_dir}/scripts"
  run "cd #{working_dir}/scripts && curl #{git_url}/scripts/csv_to_iranges.R > csv_to_iranges.R"
  run "sudo chmod +x #{working_dir}/scripts/csv_to_iranges.R"
  run "cd #{mount_point}/ESC/PET &&  #{working_dir}/scripts/csv_to_iranges.R #{working_dir}/lib/mm8ToMm9.over.chain"
  run "cd #{mount_point}/NS5/PET &&  #{working_dir}/scripts/csv_to_iranges.R #{working_dir}/lib/mm8ToMm9.over.chain"
end 
before 'pet2iranges', 'EC2:start'



desc "get the PET data"
task "fetch_pet", :roles => group_name do
  `mkdir -p results/ESC/PET`
  `mkdir -p results/NS5/PET`
  download("#{mount_point}/ESC/PET/RangedData.R", "results/ESC/PET/RangedData.R")
  download("#{mount_point}/ESC/PET/RangedData.csv", "results/ESC/PET/RangedData.csv")
  download("#{mount_point}/NS5/PET/RangedData.R", "results/NS5/PET/RangedData.R")
  download("#{mount_point}/NS5/PET/RangedData.csv", "results/NS5/PET/RangedData.csv")
end 
before 'fetch_pet', 'EC2:start'


#now we have to deal with all the  nearest stuff?

#and the tfbs stuff.

#sigh.
