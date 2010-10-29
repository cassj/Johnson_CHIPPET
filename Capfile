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
set :ebs_zone, 'eu-west-1a'  #wherever your ami is. 
set :dev, '/dev/sdf'
set :mount_point, '/mnt/data'

set :git_url, 'http://github.com/cassj/Johnson_CHIPPET/raw/master'

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


desc "Fetch expression data from Buckley Lab server"
task :expression_data, :roles => group_name do  
  run "mkdir -p #{mount_point}/publication"
  upload("data/E14-REST_DN_QC_List-w-Replicates.csv", "#{mount_point}/publication/ESC-REST-DN-raw.csv")
  upload("data/NS5-REST-DN-raw.csv", "#{mount_point}/publication/NS5-REST-DN-raw.csv")
end
before 'expression_data', 'EC2:start'




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

desc "Run limma on the differential expression data"
task :limma, :roles => group_name do
  run "mkdir -p #{working_dir}/scripts"
  
end
before 'limma', 'EC2:start'




### Process the PET data ###



