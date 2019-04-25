for i in range(45, 50):
    fout = "hs_correct_%d.sh" % i
    #print(fout)
    f_out = open(fout, "w")
    f_out.write("#$ -pe multi 16\n")
    f_out.write("#$ -cwd\n")
    #f_out.write("#$ -t 1-16 -tc 4\n")
    f_out.write("module load oraclejava/jdk1.8.0_74\n")
    f_out.write("cd /storage/home/users/pjt6/newton/newton_using_EC_reads/correction/1-overlapper\n")
    f_out.write("mkdir done_shells\n")
    f_out.write("cd correction/1-overlapper\n")
    #cmd = "./precompute.sh %d > ./overlap.0000%d.out 2>&1\n" % (i, i)
    cmd = "./precompute.sh %d > ./precompute.0000%d.out 2>&1\n" % (i, i)
    f_out.write(cmd)
    mvcm = "mv hs_correct_%d.sh ./done_shells" % i
    f_out.close()

