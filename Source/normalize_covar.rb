#!/usr/bin/env ruby

require 'narray'
FMT="%15.10f"

infile = ARGV.first

unless File.exists? infile
  $stderr.puts "ERROR: need inputfile"
  exit(1)
end

diag_file = infile + ".diag" 
ncov_file = infile + ".ncov"

t1 = []
t2 = []
ca = []  # covariance array
diag = []  # diagonal array

$stderr.puts "Reading: '#{infile}'"
File.foreach(infile) do |line| 
  c = line.strip.split
  next unless c.size == 3

  t1 << c[0]
  t2 << c[1]
  ca << c[2].to_f
  if c[0] == c[1]
    diag << Math.sqrt(c[2].to_f)
  end

end

size = Math.sqrt(t1.size).to_i

raise "Input error: diagonal terms not equal to assumed size!" unless diag.size == size

$stderr.puts "  READ MATRIX: #{size} x #{size} (#{ca.size} elements)"

cov = NArray[*ca]
cov.reshape!(size,size)
$stderr.puts "  COVAR STATS:\n\tmin = #{cov.min}\n\tmax = #{cov.max}"

# we're going to assume the covariance matrix is symmetric

$stderr.puts "Creating normalized covariance in: '#{ncov_file}'"
$stderr.puts "Creating diagonal elements in:     '#{diag_file}'"
# sleep 2 

File.open(diag_file, "w") do |fd| 
  File.open(ncov_file, "w") do |fc| 
    size.times do |i| 
      size.times do |j| 
        cov[i,j] /= (diag[i] * diag[j])

#         f.puts "%s %s #{FMT}" % [ t1[i], t2[j], cov[i,j] ] # same as input format 
        fc.print "#{FMT}" % cov[i,j]  # matrix output (part 1/2)
      end
      fc.puts # matrix output (part 2/2)
      fd.puts FMT % diag[i]
    end
  end
end

$stderr.puts "  NORMALIZED COVAR STATS:\n\tmin = #{cov.min}\n\tmax = #{cov.max}"

