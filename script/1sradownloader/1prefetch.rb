require 'fileutils'

require 'csv'

@dir= File.dirname(File.realpath(__FILE__))

def parse(input)
    data=CSV.parse(File.read("#{@dir}/#{input}"), headers: true)
    species=data.by_col[28]
    p species
    sra_code=data.by_col[0]
    folder_name=species.map {|x| x.gsub("\s","_")}
    link= data.by_col[9]
    return species, sra_code, folder_name, link
end

def remove_empty_folder
  out=File.open("#{@dir}/error.txt",'w')
  for i in 0..species.count-1
    if (Dir.entries("#{@dir}/#{folder_name[i]}") - %w{ . .. }).empty? #remove directory if no genome is downloaded
      FileUtils.rm_rf("#{@dir}/#{folder_name[i]}")
      out.puts "#{species[i]}\n"
    end
  end
end

def prefetch(code, folder)
    %x( prefetch #{code} -o #{@dir}/#{folder}/#{code}.sra >> #{@dir}/prefetch.out 2>> #{@dir}/prefetch.err )
    %x( parallel-fastq-dump --sra-id #{@dir}/#{folder}/#{code}.sra --threads 24 --split-3 --gzip --outdir #{@dir}/#{folder}/  >> #{@dir}/fastqdump.out 2>> #{@dir}/fastqdump.err)
    %x( cd #{@dir}/#{folder} && rm #{code}.sra)
end

def creating_dir_or_not(i, species_i, species_i_1, folder)
  if i == 0 || species_i != species_i_1
    FileUtils.mkdir_p("#{@dir}/#{folder}")
  end
end

def list_empty_folder(input)
  p Dir.glob("#{@dir}/#{input}").select {|f| File.directory? f}

end

def self.run(input)
  species, code, folder_name, link=parse(input)
  for i in 0..species.count-1
    creating_dir_or_not(i, species[i], species[i-1], folder_name[i])
    prefetch(code[i],folder_name[i])
  end
end


run(ARGV[0])
