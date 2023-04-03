
# Use fit hic instead of NCHG?
# Find papaer and read it.

# Set the path to the Fit-Hi-C executable
fit_hic_executable = "path/to/fitHiC.py"

# Define the input files and parameters
hicpro_contact_matrix = "path/to/hicPro_contact_matrix.txt"
fit_hic_output = "path/to/fit_hic_output.txt"
resolution = 20000  # Set the desired resolution, e.g., 20 kb

# Run Fit-Hi-C with the specified parameters
subprocess.run([
    "python",
    fit_hic_executable,
    "-f", hicpro_contact_matrix,
    "-o", fit_hic_output,
    "-r", str(resolution),
    # Add any other required parameters
])


# Chromosome parallelization for NCHG:

with open(file in input_dir) as input_files:

    chr_dict = {}
    # or list with comprehension? tuple? set? queue? stack? tree? graph? array? Think dict.
    for file in input_files:
        # find fast way to filter to chromosomes, regex? What is the best way to grep for changes in chr string, if col[0] == chr1?
        # Can also filter on inter/intra here but whY?
        # Then find way to split the file, use temp files? The input to nchg has to be files infortunately.
        # The big thing again though is concatenating the output files into one file (per resolution), before passing the whole genome file to padj.

        # Research first: Split files into temp files based on chromosome (fastest), base temp file name on chromosome + file anme (res and exp).
        # The new files can be output to a new directory, dict or filter etc to make them per chromosome, and these files are then passed to the NCHG method by just changin the input dir.
        # The NCHG method doesn't need to change, jus the input to nchg, since we need to take the input files, and then split the, by chrom, then write to new dir, and then save the output of nchg to new dir,
        # and then concatenate the output files into one file (per resolution).

        # Pass these files to the NCHG method
        # Concatenate the output files into one file (per resolution), in the correct order.

        # How does this work with the parallelization? Can we have thread and one process pool for the method? The NCHG method needs to be called with process pool.
        # Try to implement process pool for "find siginificant interactions" method.



 File "<string>", line 1, in <module>
  File "/anaconda3/envs/hic/lib/python3.8/multiprocessing/spawn.py", line 116, in spawn_main
    exitcode = _main(fd, parent_sentinel)
  File "/opt/anaconda3/envs/hic/lib/python3.8/multiprocessing/spawn.py", line 125, in _main
    prepare(preparation_data)
  File "/opt/anaconda3/envs/hic/lib/python3.8/multiprocessing/spawn.py", line 236, in prepare
    _fixup_main_from_path(data['init_main_from_path'])
  File "/opt/anaconda3/envs/hic/lib/python3.8/multiprocessing/spawn.py", line 287, in _fixup_main_from_path
    main_content = runpy.run_path(main_path,
  File "/opt/anaconda3/envs/hic/lib/python3.8/runpy.py", line 265, in run_path
Traceback (most recent call last):
Traceback (most recent call last):
  File "<string>", line 1, in <module>
  File "<string>", line 1, in <module>
  File "/Users/GBS/opt/anaconda3/envs/hic/lib/python3.8/multiprocessing/spawn.py", line 116, in spawn_main
Traceback (most recent call last):
  File "<string>", line 1, in <module>
  File "/opt/anaconda3/envs/hic/lib/python3.8/multiprocessing/spawn.py", line 116, in spawn_main
  File "/opt/anaconda3/envs/hic/lib/python3.8/multiprocessing/spawn.py", line 116, in spawn_main
    return _run_module_code(code, init_globals, run_name,
    exitcode = _main(fd, parent_sentinel)
  File "/opt/anaconda3/envs/hic/lib/python3.8/runpy.py", line 97, in _run_module_code
  File "/opt/anaconda3/envs/hic/lib/python3.8/multiprocessing/spawn.py", line 125, in _main
    exitcode = _main(fd, parent_sentinel)
  File "/opt/anaconda3/envs/hic/lib/python3.8/multiprocessing/spawn.py", line 125, in _main
    exitcode = _main(fd, parent_sentinel)
  File "/opt/anaconda3/envs/hic/lib/python3.8/multiprocessing/spawn.py", line 125, in _main
    _run_code(code, mod_globals, init_globals,
  File "/opt/anaconda3/envs/hic/lib/python3.8/runpy.py", line 87, in _run_code
    prepare(preparation_data)
    prepare(preparation_data)
  File "/opt/anaconda3/envs/hic/lib/python3.8/multiprocessing/spawn.py", line 236, in prepare
  File "/opt/anaconda3/envs/hic/lib/python3.8/multiprocessing/spawn.py", line 236, in prepare
    prepare(preparation_data)
  File "/opt/anaconda3/envs/hic/lib/python3.8/multiprocessing/spawn.py", line 236, in prepare
    exec(code, run_globals)
  File "/pymaster/script_parallel.py", line 195, in <module>
    _fixup_main_from_path(data['init_main_from_path'])
  File "/opt/anaconda3/envs/hic/lib/python3.8/multiprocessing/spawn.py", line 287, in _fixup_main_from_path
    _fixup_main_from_path(data['init_main_from_path'])
  File "/opt/anaconda3/envs/hic/lib/python3.8/multiprocessing/spawn.py", line 287, in _fixup_main_from_path
    _fixup_main_from_path(data['init_main_from_path'])
    main_content = runpy.run_path(main_path,
  File "/opt/anaconda3/envs/hic/lib/python3.8/runpy.py", line 265, in run_path
  File "/opt/anaconda3/envs/hic/lib/python3.8/multiprocessing/spawn.py", line 287, in _fixup_main_from_path
    main_content = runpy.run_path(main_path,
  File "/opt/anaconda3/envs/hic/lib/python3.8/runpy.py", line 265, in run_path
    shutil.rmtree(pbt_temp_dir)
  File "/opt/anaconda3/envs/hic/lib/python3.8/shutil.py", line 740, in rmtree
    return _run_module_code(code, init_globals, run_name,
  File "/opt/anaconda3/envs/hic/lib/python3.8/runpy.py", line 97, in _run_module_code
    main_content = runpy.run_path(main_path,
  File "/opt/anaconda3/envs/hic/lib/python3.8/runpy.py", line 265, in run_path
    return _run_module_code(code, init_globals, run_name,
  File "/opt/anaconda3/envs/hic/lib/python3.8/runpy.py", line 97, in _run_module_code
Traceback (most recent call last):
  File "<string>", line 1, in <module>
  File "/opt/anaconda3/envs/hic/lib/python3.8/multiprocessing/spawn.py", line 116, in spawn_main
    _run_code(code, mod_globals, init_globals,
  File "/opt/anaconda3/envs/hic/lib/python3.8/runpy.py", line 87, in _run_code
    return _run_module_code(code, init_globals, run_name,
  File "/opt/anaconda3/envs/hic/lib/python3.8/runpy.py", line 97, in _run_module_code
Traceback (most recent call last):
  File "<string>", line 1, in <module>
    _run_code(code, mod_globals, init_globals,
    exec(code, run_globals)
  File "/opt/anaconda3/envs/hic/lib/python3.8/runpy.py", line 87, in _run_code
  File "/pymaster/script_parallel.py", line 195, in <module>
    exitcode = _main(fd, parent_sentinel)
  File "/opt/anaconda3/envs/hic/lib/python3.8/multiprocessing/spawn.py", line 125, in _main
  File "/opt/anaconda3/envs/hic/lib/python3.8/multiprocessing/spawn.py", line 116, in spawn_main
    _run_code(code, mod_globals, init_globals,
  File "/opt/anaconda3/envs/hic/lib/python3.8/runpy.py", line 87, in _run_code
    shutil.rmtree(pbt_temp_dir)
  File "/opt/anaconda3/envs/hic/lib/python3.8/shutil.py", line 740, in rmtree
    prepare(preparation_data)
  File "/opt/anaconda3/envs/hic/lib/python3.8/multiprocessing/spawn.py", line 236, in prepare
    exitcode = _main(fd, parent_sentinel)
  File "/opt/anaconda3/envs/hic/lib/python3.8/multiprocessing/spawn.py", line 125, in _main
Traceback (most recent call last):
  File "<string>", line 1, in <module>
    exec(code, run_globals)
  File "/pymaster/script_parallel.py", line 195, in <module>
    exec(code, run_globals)
  File "/pymaster/script_parallel.py", line 195, in <module>
  File "/opt/anaconda3/envs/hic/lib/python3.8/multiprocessing/spawn.py", line 116, in spawn_main
    prepare(preparation_data)
  File "/opt/anaconda3/envs/hic/lib/python3.8/multiprocessing/spawn.py", line 236, in prepare
    exitcode = _main(fd, parent_sentinel)
  File "/opt/anaconda3/envs/hic/lib/python3.8/multiprocessing/spawn.py", line 125, in _main
    _fixup_main_from_path(data['init_main_from_path'])
  File "/opt/anaconda3/envs/hic/lib/python3.8/multiprocessing/spawn.py", line 287, in _fixup_main_from_path
    shutil.rmtree(pbt_temp_dir)
  File "/opt/anaconda3/envs/hic/lib/python3.8/shutil.py", line 740, in rmtree
    shutil.rmtree(pbt_temp_dir)
  File "/opt/anaconda3/envs/hic/lib/python3.8/shutil.py", line 740, in rmtree
    prepare(preparation_data)
  File "/opt/anaconda3/envs/hic/lib/python3.8/multiprocessing/spawn.py", line 236, in prepare
    main_content = runpy.run_path(main_path,
  File "/opt/anaconda3/envs/hic/lib/python3.8/runpy.py", line 265, in run_path
    return _rmtree_unsafe(path, onerror)
  File "/opt/anaconda3/envs/hic/lib/python3.8/shutil.py", line 618, in _rmtree_unsafe
    _fixup_main_from_path(data['init_main_from_path'])
  File "/opt/anaconda3/envs/hic/lib/python3.8/multiprocessing/spawn.py", line 287, in _fixup_main_from_path
    return _run_module_code(code, init_globals, run_name,
  File "/anaconda3/envs/hic/lib/python3.8/runpy.py", line 97, in _run_module_code
    return _rmtree_unsafe(path, onerror)
  File "/opt/anaconda3/envs/hic/lib/python3.8/shutil.py", line 618, in _rmtree_unsafe
    return _rmtree_unsafe(path, onerror)
  File "/opt/anaconda3/envs/hic/lib/python3.8/shutil.py", line 618, in _rmtree_unsafe
    main_content = runpy.run_path(main_path,
  File "/opt/anaconda3/envs/hic/lib/python3.8/runpy.py", line 265, in run_path
    _run_code(code, mod_globals, init_globals,
    onerror(os.unlink, fullname, sys.exc_info())
  File "/opt/anaconda3/envs/hic/lib/python3.8/runpy.py", line 87, in _run_code
  File "/opt/anaconda3/envs/hic/lib/python3.8/shutil.py", line 616, in _rmtree_unsafe
    exec(code, run_globals)
  File "/Users/GBS/Master/pymaster/script_parallel.py", line 195, in <module>
    return _run_module_code(code, init_globals, run_name,
  File "/opt/anaconda3/envs/hic/lib/python3.8/runpy.py", line 97, in _run_module_code
    onerror(os.unlink, fullname, sys.exc_info())
  File "/opt/anaconda3/envs/hic/lib/python3.8/shutil.py", line 616, in _rmtree_unsafe
    onerror(os.unlink, fullname, sys.exc_info())
  File "/opt/anaconda3/envs/hic/lib/python3.8/shutil.py", line 616, in _rmtree_unsafe
    shutil.rmtree(pbt_temp_dir)
    _run_code(code, mod_globals, init_globals,
  File "/opt/anaconda3/envs/hic/lib/python3.8/runpy.py", line 87, in _run_code
    os.unlink(fullname)
    return _rmtree_unsafe(path, onerror)
FileNotFoundError: [Errno 2] No such file or directory: '/Pipeline_out/chr18_INC/chr18_norm/temp_dir/pybedtools_tmp/pybedtools.estq_cc6.tmp'
  File "/opt/anaconda3/envs/hic/lib/python3.8/shutil.py", line 618, in _rmtree_unsafe
  File "/opt/anaconda3/envs/hic/lib/python3.8/shutil.py", line 740, in rmtree
    exec(code, run_globals)
  File "/pymaster/script_parallel.py", line 195, in <module>
    os.unlink(fullname)
FileNotFoundError: [Errno 2] No such file or directory: '/Pipeline_out/chr18_INC/chr18_norm/temp_dir/pybedtools_tmp/pybedtools.01qt0pk6.tmp'
    os.unlink(fullname)
FileNotFoundError: [Errno 2] No such file or directory: '/Pipeline_out/chr18_INC/chr18_norm/temp_dir/pybedtools_tmp/pybedtools.estq_cc6.tmp'
    shutil.rmtree(pbt_temp_dir)
  File "/opt/anaconda3/envs/hic/lib/python3.8/shutil.py", line 740, in rmtree
    return _rmtree_unsafe(path, onerror)
  File "/opt/anaconda3/envs/hic/lib/python3.8/shutil.py", line 618, in _rmtree_unsafe
    onerror(os.unlink, fullname, sys.exc_info())
  File "/opt/anaconda3/envs/hic/lib/python3.8/shutil.py", line 616, in _rmtree_unsafe
    return _rmtree_unsafe(path, onerror)
  File "/opt/anaconda3/envs/hic/lib/python3.8/shutil.py", line 618, in _rmtree_unsafe
    onerror(os.unlink, fullname, sys.exc_info())
  File "/opt/anaconda3/envs/hic/lib/python3.8/shutil.py", line 616, in _rmtree_unsafe
    onerror(os.unlink, fullname, sys.exc_info())
  File "/opt/anaconda3/envs/hic/lib/python3.8/shutil.py", line 616, in _rmtree_unsafe
    os.unlink(fullname)
    os.unlink(fullname)
FileNotFoundError: [Errno 2] No such file or directory: '/HiC-Data/Pipeline_out/chr18_INC/chr18_norm/temp_dir/pybedtools_tmp/pybedtools.8x3nykj_.tmp'
FileNotFoundError: [Errno 2] No such file or directory: '/HiC-Data/Pipeline_out/chr18_INC/chr18_norm/temp_dir/pybedtools_tmp/pybedtools.bz2si1z8.tmp'
    os.unlink(fullname)
FileNotFoundError: [Errno 2] No such file or directory: '/HiC-Data/Pipeline_out/chr18_INC/chr18_norm/temp_dir/pybedtools_tmp/pybedtools.xj865afe.tmp'
    _fixup_main_from_path(data['init_main_from_path'])
  File "/opt/anaconda3/envs/hic/lib/python3.8/multiprocessing/spawn.py", line 287, in _fixup_main_from_path
    main_content = runpy.run_path(main_path,
  File "/opt/anaconda3/envs/hic/lib/python3.8/runpy.py", line 265, in run_path
    return _run_module_code(code, init_globals, run_name,
  File "/opt/anaconda3/envs/hic/lib/python3.8/runpy.py", line 97, in _run_module_code
    _run_code(code, mod_globals, init_globals,
  File "/opt/anaconda3/envs/hic/lib/python3.8/runpy.py", line 87, in _run_code
    exec(code, run_globals)
  File "/pymaster/script_parallel.py", line 195, in <module>
    shutil.rmtree(pbt_temp_dir)
  File "/opt/anaconda3/envs/hic/lib/python3.8/shutil.py", line 740, in rmtree
    return _rmtree_unsafe(path, onerror)
  File "/opt/anaconda3/envs/hic/lib/python3.8/shutil.py", line 618, in _rmtree_unsafe
    onerror(os.unlink, fullname, sys.exc_info())
  File "/opt/anaconda3/envs/hic/lib/python3.8/shutil.py", line 616, in _rmtree_unsafe
    os.unlink(fullname)
FileNotFoundError: [Errno 2] No such file or directory: '/HiC-Data/Pipeline_out/chr18_INC/chr18_norm/temp_dir/pybedtools_tmp/pybedtools.up48ly8a.tmp'
 43%|█████████████████████████████████████████████████████████▍                                                                            | 3/7 [00:21<00:29,  7.26s/it]
Traceback (most recent call last):
  File "script_parallel.py", line 848, in <module>
    run_pipeline()
  File "script_parallel.py", line 841, in run_pipeline
    method()  # Call the method
  File "script_parallel.py", line 610, in input_to_nchg
    futures = list(executor.map(Pipeline.find_siginificant_interactions, full_paths))
  File "/opt/anaconda3/envs/hic/lib/python3.8/concurrent/futures/process.py", line 484, in _chain_from_iterable_of_lists
    for element in iterable:
  File "/opt/anaconda3/envs/hic/lib/python3.8/concurrent/futures/_base.py", line 619, in result_iterator
    yield fs.pop().result()
  File "/opt/anaconda3/envs/hic/lib/python3.8/concurrent/futures/_base.py", line 444, in result
    return self.__get_result()
  File "/opt/anaconda3/envs/hic/lib/python3.8/concurrent/futures/_base.py", line 389, in __get_result
    raise self._exception
concurrent.futures.process.BrokenProcessPool: A process in the process pool was terminated abruptly while the future was running or pending.


