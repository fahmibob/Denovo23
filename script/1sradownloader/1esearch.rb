require 'csv'
require 'fileutils'

class Wgs
  def self.run(input) #it takes csv file as input
    dir= File.dirname(File.realpath(__FILE__)) #get directory of the file from which this method is called
    data=CSV.parse(File.read("#{dir}/#{input}"), headers: true)
    query=data.by_col[7] #get the ftp url in the second column of defined csv file

    for i in 0..query.count-1
      %x( sh ./1esearch.sh #{query[i]}) #download genome and save to each directory
    end
  end
end

Wgs.run(ARGV[0])
