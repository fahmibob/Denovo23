def checkDuplicate(glob_1,glob_2)
  dir1=Dir.glob(glob_1)
  dir2=Dir.glob(glob_2)
  dir1=dir1.map{|x| x.split("/")[-1]}
  dir2=dir2.map{|x| x.split("/")[-1]}
  k=dir1.&(dir2)
  puts k
end
