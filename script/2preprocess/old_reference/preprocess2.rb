require 'csv'
require 'fileutils'
require 'pathname'
require 'zip'
require 'json'
require_relative 'module.rb'

$dir=File.expand_path("..", Dir.pwd)
$dirSpeciesDef="#{$dir}/species"
$logdir="#{$dir}/log"
$outdir="#{$dir}/output"
$database=CSV.parse(File.read("#{$dir}/data/run_list.csv"), headers: true)

class Input
  include Dir_utils, Data_utils
  attr_reader :inputDir, :listSpecies
  def initialize(inputDir, listSpecies=nil)
    @inputDir=inputDir
    @listSpecies=listSpecies
  end

  def single_double_mode(filteredSingleEnds, filteredSingleEndsName, pairEnds, pairEndsName)
    single_or_double_end=[]
    files_location=[]
    files_name=[]
    filteredSingleEnds.each_with_index do |indv, index|
      single_or_double_end.push(0)
      files_location.push(indv)
      files_name.push(filteredSingleEndsName[index])
    end

    pairEnds.each_with_index do |indv, index|
      single_or_double_end.push(1)
      files_location.push(indv)
      files_name.push(pairEndsName[index])
    end
    return files_name, files_location, single_or_double_end
  end

  def whole_data(filteredSingleEnds, filteredSingleEndsName, pairEnds, pairEndsName)
    files_location=[]
    files_name=[]
    filteredSingleEnds.each_with_index do |indv, index|
      files_location.push(indv)
      files_name.push(filteredSingleEndsName[index])
    end
    pairEnds.each_with_index do |indv, index|
      files_name.push(pairEndsName[index])
      files_location.push(pairForward(pairEnds[index]))
      files_name.push(pairEndsName[index])
      files_location.push(indv)
    end
    return files_name, files_location
  end

  def parse(separate_forward_reverse, filetype) #mode 0 separate forward reverse analysis (fastp); 1 join both type of analysis (fastqc)
    if Dir.exist?(Dir.glob("#{@inputDir}/*/*")[0])
      inputFiles=collect_data("#{@inputDir}/*", "*#{filetype}").flatten
      platform=true #platformwise
    else
      inputFiles=collect_data(@inputDir, "*#{filetype}").flatten
      platform=false #non_platformwise
    end

    unless @listSpecies.nil?
      speciesQuery=[]
      selectedfiles=[]
      File.readlines(@listSpecies).each {|line|
        species=line.gsub("\n","")
        species=species.gsub("\s","_")
        speciesQuery.push(species)
      }
      speciesOfInputFiles=inputFiles.map {|file| detect_species(file)}

      for k in  0..speciesQuery.count-1
        for i in 0..speciesOfInputFiles.count-1
          if speciesOfInputFiles[i] == speciesQuery[k]
            selectedfiles.push(inputFiles[i])
          end
        end
      end
      inputFiles=selectedfiles
    end


    singleEnds=inputFiles.select {|file| file.match(/\d+\d+.fastq.gz/)}
    pairEnds=inputFiles.select {|file| file.match("_2.fastq.gz")}
    pairEndsName=pairEnds.map {|file| file.split("/").last.match("_2.fastq.gz").pre_match}
    singleEndsName=singleEnds.map {|file| file.split("/").last.match(".fastq.gz").pre_match}
    filteredSingleEndsName=singleEndsName.-(pairEndsName)
    filteredSingleEnds=filteredSingleEndsName.map {|indv| singleEnds[singleEndsName.index(indv)]}

    if separate_forward_reverse
      result=single_double_mode(filteredSingleEnds, filteredSingleEndsName, pairEnds, pairEndsName)
      files_name, files_location=result[0], result[1]
      single_or_double_end=result[2]
      return files_name, files_location, single_or_double_end, platform
    else
      result=whole_data(filteredSingleEnds, filteredSingleEndsName, pairEnds, pairEndsName)
      files_name, files_location=result[0], result[1]
      return files_name, files_location, platform
    end
  end
end

class Log
  def self.time_diff(start_time, end_time)
    seconds_diff = (start_time - end_time).to_i.abs

    hours = seconds_diff / 3600
    seconds_diff -= hours * 3600

    minutes = seconds_diff / 60
    seconds_diff -= minutes * 60

    seconds = seconds_diff

    "#{hours.to_s.rjust(2, '0')}:#{minutes.to_s.rjust(2, '0')}:#{seconds.to_s.rjust(2, '0')}"

  end

  def self.file_size(file)
    return File.size("#{file}").to_f / 2**20
  end

  def self.libraryType_size(librarytype, file)
    if librarytype == "pair"
      return (file_size(file)*2).round(2)
    elsif librarytype == "single"
      return file_size(file).round(2)
    end
  end

  def self.timelog_header(outputlog)
    header="Species\tCode\tsize(MB)\tProcessing_time"
    if File.file?(outputlog)
      unless File.readlines(outputlog)[0].match("Species\tCode")
        write=File.open(outputlog, "a+")
        write.puts header
        write.close
      end
    else
      write=File.open(outputlog, "a+")
      write.puts header
      write.close
    end
  end

  def self.timelog(starting, ending, librarytype, filedir, filename, outputlog) #type paired = true, single = false
    elapsed = time_diff(starting, ending)
    if elapsed != "00:00:00"
      size=libraryType_size(librarytype, filedir)
      note="#{filedir.split("/")[-2]}\t#{filename}\t#{size}\t#{elapsed}"
      write=File.open(outputlog, "a+")
      write.puts note
      write.close
    end
  end
end

class Fastp
  include ReadDatabase, Dir_utils, Data_utils
  attr_reader :jobName, :inputDir, :listSpecies

  def initialize(jobName, inputDir, listSpecies=nil)
    @jobName=jobName
    @logfile="#{$logdir}/#{@jobName}_fastp.log"
    @outdir="#{$outdir}/fastp"
    @listSpecies=listSpecies
    if inputDir.nil?
      @inputDir=$dirSpeciesDef
    else
      @inputDir=inputDir
    end
  end

  def execute(t=10, platformwiseOut=false)
    outputdir="#{@outdir}/#{@jobName}"
    Dir.mkdir(outputdir) unless Dir.exist?(outputdir)

    inputData=Input.new(@inputDir, @listSpecies).parse(true, "fastq.gz")

    fileName, inputLocation =inputData[0], inputData[1]
    singleOrPairEnds, platformwiseIn =inputData[2], inputData[3]

    fileName.each_with_index do |file, index|
      platform=check_platform(platformwiseIn, platformwiseOut, inputLocation[index], file, outputdir)
      species=detect_species(inputLocation[index])
      Dir.mkdir("#{outputdir}/#{platform}/#{species}") \
      unless Dir.exist?("#{outputdir}/#{platform}/#{species}")
      html="#{outputdir}/report/#{platform}/#{file}_fastp.html"
      json="#{outputdir}/report/#{platform}/#{file}_fastp.json"
      Dir.mkdir("#{outputdir}/report") unless Dir.exist?("#{outputdir}/report")
      Dir.mkdir("#{outputdir}/report/#{platform}") \
      unless Dir.exist?("#{outputdir}/report/#{platform}")

      if singleOrPairEnds[index] == 0 #single end
        output="#{outputdir}/#{platform}/#{species}/FP#{file}.fastq.gz"
        puts "processing #{file} for fastp"
        unless Dir.glob("#{output}").any?
          unless system("bash ./tools.sh fastp_single \
            #{inputLocation[index]} #{output} #{html} \
            #{json} #{t}")
            puts "There was an error for #{file}"
          end
        end
      else #paired ends
        input_f=pairForward(inputLocation[index])
        input_r=pairReverse(inputLocation[index])
        output_f="#{outputdir}/#{platform}/#{species}/FP#{file}_1.fastq.gz"
        output_r="#{outputdir}/#{platform}/#{species}/FP#{file}_2.fastq.gz"
        unless Dir.glob("#{output_f}").any?
          puts "processing #{file} for fastp"
          unless system("bash ./tools.sh fastp_paired #{input_f} \
            #{input_r} #{output_f} #{output_r} #{html} #{json} #{t}")
            puts "There was an error for #{name}"
          end
        end
      end
    end
  end
end

class Fastqc
  include ReadDatabase, Dir_utils, Data_utils
  attr_reader :jobName, :inputDir, :listSpecies

  def initialize(jobName, inputDir, listSpecies=nil)
    @jobName=jobName
    @logfile="#{$logdir}/#{@jobName}_fastqc.log"
    @outdir="#{$outdir}/fastqc"
    @listSpecies=listSpecies
    if inputDir.nil?
      @inputDir=$dirSpeciesDef
    else
      @inputDir=inputDir
    end
  end

  def fastqc_join_array(files, joined_number=5)
    if files.count==0 || files.count==1
      input=files
    elsif files.count/joined_number == 0
      input=files.join(" ")
    else
      fives=files.count/joined_number
      input=[]

      for i in 1..fives
        input.push(files[(i-1)*joined_number..((i-1)*joined_number)+(joined_number-1)].join(" "))
      end

      if files.count%joined_number != 0
        input.push(files[-(files.count%joined_number)..-1].join(" "))
      end
    end
    return input
  end

  def fastqc_platformwise(array, t, log, outputdir)
    fastqcInputs=fastqc_join_array(array.drop(1),5)
    fastqcInputs.each do |input|
      puts "processing #{input}"
      unless system("fastqc #{input} --outdir \
        #{outputdir}/#{array[0]} -threads #{t} 2>>#{@logfile}")
        puts "There was an error for #{input}"
      end
    end
  end

  def execute(t=10, platformwiseOut=false)
    outputdir="#{@outdir}/#{@jobName}"
    Dir.mkdir(outputdir) unless Dir.exist?(outputdir)

    inputData=Input.new(@inputDir, @listSpecies).parse(false, "fastq.gz")

    fileName=inputData[0]
    inputLocation=inputData[1]
    platformwiseIn=inputData[2]

    platforms=[]
    illumina, bgi=["ILLUMINA"],["BGISEQ"]
    abi_solid, ion, ls454=["ABI_SOLID"],["ION_TORRENT"],["LS454"]
    fileName.each_with_index { |file, index|
      platform=check_platform(platformwiseIn, platformwiseOut, inputLocation[index], file, outputdir)
      platforms.push(platform)
      if platform=="ILLUMINA"
        illumina.push(inputLocation[index])
      elsif platform=="LS454"
        ls454.push(inputLocation[index])
      elsif platform=="ABI_SOLID"
        abi_solid.push(inputLocation[index])
      elsif platform=="BGISEQ"
        bgi.push(inputLocation[index])
      elsif platform=="ION_TORRENT"
        ion.push(inputLocation[index])
      end
    }

    if platformwiseIn || platformwiseOut
      for i in 0..platforms.uniq.count-1
        if platforms.uniq[i]==illumina[0]
          fastqc_platformwise(illumina, t, @logfile, outputdir)
        elsif platforms.uniq[i]==ls454[0]
          fastqc_platformwise(ls454, t, @logfile, outputdir)
        elsif platforms.uniq[i]==abi_solid[0]
          fastqc_platformwise(abi_solid, t, @logfile, outputdir)
        elsif platforms.uniq[i]==bgi[0]
          fastqc_platformwise(bgi, t, @logfile, outputdir)
        elsif platforms.uniq[i]==ion[0]
          fastqc_platformwise(ion, t, @logfile, outputdir)
        end
      end
    else
      fastqcInputs=fastqc_join_array(inputLocation, 5)
      fastqcInputs.each do |input|
        unless system("fastqc #{input} --outdir \
          #{outputdir} -threads #{t} 2>>#{@logfile}")
          puts "There was an error for #{input}"
        end
      end
    end
  end


end

class Multiqc
  include Dir_utils
  attr_reader :directory

  def initialize(directory)
    @directory=directory
    @multiqc="#{$dir}/output/multiqc"
  end

  def dividing_multiqc(files_amount)
    maximum_capacity=500
    fraction=files_amount / maximum_capacity
    index_i, index_n=[],[]
    if files_amount <= 500
      index_i.push(0)
      index_n.push(files_amount-1)
    elsif files_amount % maximum_capacity == 0 && files_amount > 500
      first_index=0
      for i in 1..fraction
        last_index=first_index+maximum_capacity-1
        index_i.push(first_index)
        index_n.push(last_index)
        first_index+=maximum_capacity
      end
    elsif files_amount % maximum_capacity != 0
      first_index=0
      opt_frac=files_amount/(files_amount/maximum_capacity.to_f).ceil
      for i in 1..fraction
        last_index=first_index+opt_frac-1
        index_i.push(first_index)
        index_n.push(last_index)
        first_index+=opt_frac
      end
      rest=files_amount-last_index
      index_i.push(last_index+1)
      index_n.push(files_amount-1)
    end
    return index_i,index_n
  end

  def sequence_length_fastp(json)
    file=File.read(json)
    data_hash=JSON.parse(file)
    length=data_hash["summary"]["before_filtering"]["read1_mean_length"]
    return length.to_i
  end

  def sequence_length(zipfile_dir)
    seq=String.new
    Zip::File.open(zipfile_dir) do |files|
      files.each do |entry|
        if entry.to_s.match("/").post_match == "fastqc_data.txt"
          content=entry.get_input_stream.read
          seqlength=content.match("Sequence length\t").post_match.match("\n").pre_match
          seq=seqlength
          if seqlength.match("-")
            seqlength=seqlength.match("-").post_match
            seq=seqlength
          end
        end
      end
    end
    return seq.to_i
  end

  def running_multiqc(data, type)
    if data.all?
      ary=dividing_multiqc(data.count)
      dir_name=@directory.split("/").last
      output="#{@multiqc}/#{dir_name}"
      Dir.mkdir(output) unless File.exists?(output)
      for i in 0..ary[0].count-1
        File.delete("#{output}/#{type}#{i+1}.txt") if File.exist?("#{output}/#{type}#{i+1}.txt")
        write=File.open("#{output}/#{type}#{i+1}.txt","a+")
        for k in ary[0][i]..ary[1][i]
          write.puts data[k]
        end
        write.close
      end

      for i in 0..ary[0].count-1
        unless system("bash ./tools.sh multi #{output}/#{type}#{i+1}.txt")
          puts "There was an error for #{output}/#{type}#{i+1}.txt"
        end
      end

      results=Dir.glob("#{File.dirname(File.realpath(__FILE__))}/multiqc*")
      results.each {|f| FileUtils.mv(f, "#{output}/#{f.split('/').last}")}
    else
      puts "#{data} has no content"
    end

  end

  def multiqc_fastp
    data=collect_data_in_folder(@directory,"*.json")
    short,mid_low,mid_high,long=[],[],[],[]
    data.each {|f|
      length=sequence_length_fastp(f)
      if length <150
        short.push(f)
      elsif length >=150 && length <750
        mid_low.push(f)
      elsif length >=750 && length <3750
        mid_high.push(f)
      else
        long.push(f)
      end
    }

    running_multiqc(short,"short") if short.any?
    running_multiqc(mid_low,"mid_low") if mid_low.any?
    running_multiqc(mid_high,"mid_high") if mid_high.any?
    running_multiqc(long,"long") if long.any?

  end

  def multiqc_fastqc
    data=collect_data_in_folder(@directory,"*.zip")

    short,mid_low,mid_high,long=[],[],[],[]
    data.each {|f|
      if sequence_length(f) <150
        short.push(f)
      elsif sequence_length(f) >=150 && sequence_length(f) <750
        mid_low.push(f)
      elsif sequence_length(f) >=750 && sequence_length(f) <3750
        mid_high.push(f)
      else
        long.push(f)
      end
    }

    running_multiqc(short,"short") if short.any?
    running_multiqc(mid_low,"mid_low") if mid_low.any?
    running_multiqc(mid_high,"mid_high") if mid_high.any?
    running_multiqc(long,"long") if long.any?

  end
end

class Sortmerna
  include ReadDatabase, Dir_utils, Data_utils
  attr_reader :jobName, :inputDir, :listSpecies

  def initialize(jobName, inputDir, listSpecies=nil)
    @jobName=jobName
    @outdir="#{$outdir}/sortmerna"
    @kvdb="/home/fahmi/sortmerna/run/kvdb"
    @listSpecies=listSpecies
    if inputDir.nil?
      @inputDir=$dirSpeciesDef
    else
      @inputDir=inputDir
    end
  end

  def execute(t=10, parallel=nil, platformwiseOut=false)
    sleep(10)

    outputdir="#{@outdir}/#{@jobName}"
    logfile="#{$logdir}/#{@jobName}_sortme#{parallel}.log"
    timelog="#{$logdir}/#{@jobName}_time_sortme#{parallel}.log"

    Dir.mkdir(outputdir) unless Dir.exist?(outputdir)

    inputData=Input.new(@inputDir, @listSpecies).parse(true, ".fastq.gz")
    fileName, inputLocation =inputData[0], inputData[1]
    singleOrPairEnds, platformwiseIn =inputData[2], inputData[3]

    Log.timelog_header(timelog)

    fileName.each_with_index do |file, index|
      platform=check_platform(platformwiseIn, platformwiseOut, inputLocation[index], file, outputdir)
      species=detect_species(inputLocation[index])
      Dir.mkdir("#{outputdir}/#{platform}/#{species}") \
      unless Dir.exist?("#{outputdir}/#{platform}/#{species}")

      if singleOrPairEnds[index] == 0 #single end
        starting=Time.now
        output="#{outputdir}/#{platform}/#{species}/SMR#{file}"
        unless Dir.glob("#{output}.fq.gz").any?
          puts "processing #{species} for sortmerna"
          unless system("bash ./tools.sh singleSortRna \
            #{inputLocation[index]} #{output} #{t} #{parallel} 1>> #{logfile}" )
            puts "There was an error for #{species} #{file}"
          end
        end
      ending=Time.now
        Log.timelog(starting, ending, "single", inputLocation[index], file, timelog)
      end

      if singleOrPairEnds[index] == 1 #paired end
        starting=Time.now
        input_f=pairForward(inputLocation[index])
        input_r=pairReverse(inputLocation[index])
        output="#{outputdir}/#{platform}/#{species}/SMR#{file}"
        unless Dir.glob("#{output}_fwd.fq.gz").any?
          puts "processing #{species} for sortmerna"
          unless system("bash ./tools.sh pairSortRna #{input_f} \
            #{input_r} #{output} #{t} #{parallel} 1>> #{logfile}")
            puts "There was an error for #{species} #{file}"
          end
        end
        ending=Time.now
        Log.timelog(starting, ending, "pair", input_f, file, timelog)
      end

    end
  end
end
