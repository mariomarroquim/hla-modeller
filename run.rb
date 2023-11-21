require 'concurrent'

allele_sequences_file = 'hla_prot.fasta'
current_allele_name = nil
allele_sequences = {}
skip_line = false

ALLOWED_LOCI = %w[A B C DRB1 DRB3 DRB4 DRB5 DQA1 DQB1 DPA1 DPB1]
ESM_TIMEOUT = 32

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

signal_peptide_sizes = {
  'A' => 24,
  'B' => 24,
  'C' => 24,
  'DRB1' => 29,
  'DRB3' => 29,
  'DRB4' => 29,
  'DRB5' => 29,
  'DQA1' => 23,
  'DQB1' => 32,
  'DPA1' => 31,
  'DPB1' => 29
}

mature_sequence_sizes = {
  'A' => 274,
  'B' => 276,
  'C' => 274,
  'DRB1' => 190,
  'DRB3' => 190,
  'DRB4' => 190,
  'DRB5' => 190,
  'DQA1' => 183,
  'DQB1' => 189,
  'DPA1' => 181,
  'DPB1' => 189
}

filtered_allele_sequences = {}

allele_sequences.each do |allele, sequence|
  locus = allele.split('*').first

  if sequence.size < mature_sequence_sizes[locus]
    next

  elsif (sequence.size + signal_peptide_sizes[locus]) >= mature_sequence_sizes[locus]
    mature_sequence = sequence[signal_peptide_sizes[locus]...(mature_sequence_sizes[locus] + signal_peptide_sizes[locus])]

  else
    mature_sequence = sequence[0...(mature_sequence_sizes[locus] - 1)]
  end

  filtered_allele_sequences[allele] = mature_sequence
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
  pool.post do
    `mkdir -p output/HLA_#{allele.split('*').first}`

    `curl -X POST -s --connect-timeout #{ESM_TIMEOUT} --data "#{sequence}" https://api.esmatlas.com/foldSequence/v1/pdb/ > output/HLA_#{allele.split('*').first}/#{allele.gsub('*', '_').gsub(
      ':', '_'
    )}.pdb`
    sleep ESM_TIMEOUT
  end
end

pool.shutdown
pool.wait_for_termination
