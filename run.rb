require 'concurrent'

allele_sequences_file = 'hla_prot.fasta'
current_allele_name = nil
allele_sequences = {}
skip_line = false

ALLOWED_LOCI = %w[DRB1 DRB3 DRB4 DRB5 DQA1 DQB1 DPA1 DPB1]
ESM_TIMEOUT = 16

unless File.exist?(allele_sequences_file)
  system('wget', "https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/#{allele_sequences_file}")
end

File.open(allele_sequences_file).each_line do |line|
  if line.include?('>')
    current_allele_name = line.split(' ')[1].split(':')[0..1].join(':')
    current_locus = current_allele_name.split('*').first

    unless allele_sequences[current_allele_name].nil?
      skip_line = true
      next
    end

    unless ALLOWED_LOCI.include?(current_locus)
      skip_line = true
      next
    end

    if current_allele_name[-1].match(/[A-Za-z]/)
      skip_line = true
      next
    end

    allele_sequences[current_allele_name] = ''
    skip_line = false

  elsif skip_line == false
    allele_sequences[current_allele_name] += line.strip
  end
end

File.delete(allele_sequences_file)

puts "DOWNLOADED ALLELES QTD: #{allele_sequences.size}"

min_sequence_sizes = {
  'A' => 270,
  'B' => 270,
  'C' => 270,
  'DRB1' => 180,
  'DRB3' => 180,
  'DRB4' => 180,
  'DRB5' => 180,
  'DQA1' => 180,
  'DQB1' => 180,
  'DPA1' => 180,
  'DPB1' => 180
}

filtered_allele_sequences = {}

allele_sequences.each do |allele, sequence|
  locus = allele.split('*').first

  next if sequence.size < min_sequence_sizes[locus]

  filtered_allele_sequences[allele] = sequence
end

puts "FILTERED ALLELES QTD: #{filtered_allele_sequences.size}"

output = File.new('input/sequences.fasta', 'w')
output.sync = true

filtered_allele_sequences.each do |allele, sequence|
  output.puts ">#{allele}".strip

  sequence.split('').each_slice(60) do |slice|
    output.puts slice.join.strip
  end
end

output.flush
output.close

puts "MODELING WHAT IS POSSIBLE USING #{Concurrent.processor_count} THREADS..."

pool = Concurrent::FixedThreadPool.new(Concurrent.processor_count)

filtered_allele_sequences.each do |allele, sequence|
  output_folder = "output/HLA_#{allele.split('*').first}"
  model_name = "#{allele.gsub('*', '_').gsub(':', '_')}.pdb"

  next if File.exist?("#{output_folder}/#{model_name}") && File.size("#{output_folder}/#{model_name}") > 0
  `rm -f #{output_folder}/#{model_name}` if File.exist?("#{output_folder}/#{model_name}")
  
  pool.post do
    `mkdir -p #{output_folder}`

    puts "[#{Time.now}] Running: curl -X POST -s --insecure --connect-timeout #{ESM_TIMEOUT} --data \"#{sequence}\" https://api.esmatlas.com/foldSequence/v1/pdb/ > #{output_folder}/#{model_name}"

    `curl -X POST -s --insecure --connect-timeout #{ESM_TIMEOUT} --data "#{sequence}" https://api.esmatlas.com/foldSequence/v1/pdb/ > #{output_folder}/#{model_name}`

    sleep ESM_TIMEOUT
  end
end

pool.shutdown
pool.wait_for_termination
