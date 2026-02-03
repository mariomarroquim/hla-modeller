# frozen_string_literal: true

# Main settings for this script:
TARGET_LOCI = %w[A B C DRB1 DRB3 DRB4 DRB5 DQB1 DPB1].freeze
TARGET_ALLELES = File.readlines('target_alleles.txt', chomp: true)

puts 'STARTED'

allele_sequences_file = 'hla_prot.fasta'
current_allele_name = nil
allele_sequences = {}
skip_line = false

File.delete(allele_sequences_file) if File.exist?(allele_sequences_file)
url = "https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/#{allele_sequences_file}"
system('wget', '-q', url, out: File::NULL, err: File::NULL)

File.open(allele_sequences_file).each_line do |line|
  if line.include?('>')
    current_allele_name = line.split(' ')[1].split(':')[0..1].join(':')
    current_locus = current_allele_name.split('*').first

    unless allele_sequences[current_allele_name].nil?
      skip_line = true
      next
    end

    if current_allele_name[-1].match(/[A-Z]/)
      skip_line = true
      next
    end

    if !TARGET_LOCI.empty? && !TARGET_LOCI.include?(current_locus)
      skip_line = true
      next
    end

    if !TARGET_ALLELES.empty? && !TARGET_ALLELES.include?(current_allele_name.tr('*', ''))
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

puts "DOWNLOADED #{allele_sequences.size} ALLELES"

signal_peptides_size = {
  'A' => 24,
  'B' => 24,
  'C' => 24,
  'DRB1' => 29,
  'DRB3' => 29,
  'DRB4' => 29,
  'DRB5' => 29,
  'DQB1' => 32,
  'DPB1' => 29
}

sequences_size = {
  'A' => 270,
  'B' => 270,
  'C' => 270,
  'DRB1' => 180,
  'DRB3' => 180,
  'DRB4' => 180,
  'DRB5' => 180,
  'DQB1' => 180,
  'DPB1' => 180
}

filtered_allele_sequences = {}

allele_sequences.each do |allele, sequence|
  locus = allele.split('*').first

  locus_sequence_size = sequences_size[locus]
  next if sequence.size < locus_sequence_size

  locus_peptide_size = signal_peptides_size[locus]
  mature_sequence = sequence[locus_peptide_size..(locus_peptide_size + locus_sequence_size) - 1]
  next if mature_sequence.size < locus_sequence_size

  filtered_allele_sequences[allele] = mature_sequence
end

puts "SELECTED #{filtered_allele_sequences.size} ALLELES"

output = File.new('input/sequences.fasta', 'w')
output.sync = true

filtered_allele_sequences.each do |allele, sequence|
  output.puts ">#{allele}"

  sequence.split('').each_slice(60) do |slice|
    output.puts slice.join
  end
end

output.close

puts 'MODELING...'

esm_min_file_size = 41
esm_timeout = 32

filtered_allele_sequences.each do |allele, sequence|
  locus_path = "output/HLA_#{allele.split('*').first}"
  allele_path = "#{locus_path}/#{allele.gsub('*', '_').gsub(':', '_')}.pdb"

  allele_file_size = File.size?(allele_path)
  next if allele_file_size && allele_file_size > esm_min_file_size

  `mkdir -p #{locus_path}`

  `rm -f #{allele_path}` if allele_file_size

  url = 'https://api.esmatlas.com/foldSequence/v1/pdb/'

  `curl -X POST -s --insecure --connect-timeout #{esm_timeout} --data "#{sequence}" #{url} > #{allele_path}`
end

puts 'FINISHED!'
