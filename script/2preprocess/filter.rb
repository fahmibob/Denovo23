input1="l3.txt"
input2="list1.txt"
thisdir=File.dirname(__FILE__)
ref=File.readlines("#{thisdir}/#{input2}").map{|x|x.gsub("\n","").gsub("SMR","")}
query=File.readlines("#{thisdir}/#{input1}").map{|x|x.gsub("\n","")}
out=query.-(ref)

puts out.count
hehe=File.open("ls3.txt","w")
hehe.puts(out)
