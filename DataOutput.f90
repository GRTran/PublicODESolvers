module dataOutput

implicit none

contains

  subroutine odeDataWriter(input_array, nvars, ref_num, file_name)
    REAL(KIND=8), DIMENSION(:,:)    ::    input_array
    INTEGER                         ::    nvars
    INTEGER                         ::    ref_num
    CHARACTER(LEN=*)                ::    file_name
    INTEGER                         ::    ios
    INTEGER                         ::    row_counter

    if(nvars<1) stop "Error on input to odeDataWriter()"

    open(unit=ref_num, file=file_name, iostat=ios, action="write")
    if ( ios /= 0 ) stop "Error opening file "

    do row_counter=1,size(input_array(:,1))

      !! open file and write
      write(unit=ref_num, fmt=*, iostat=ios) input_array(row_counter,:)
      if ( ios /= 0 ) stop "Write error in file unit "
    enddo

    close(ref_num)
  end subroutine

  subroutine odeGnuplotDataPlotter(input_file, file_ref, output_file)
    CHARACTER(LEN=*)    ::    input_file
    INTEGER             ::    file_ref
    CHARACTER(LEN=*)    ::    output_file
    LOGICAL             ::    file_exists
    INTEGER             ::    stat


    inquire(file=input_file,exist=file_exists)
    if (file_exists) then
      open(file_ref, file = output_file)
      write(file_ref,*) "set term png"
      write(file_ref,*) "set output 'output/bench/bench_dat.png'"
      write(file_ref,*) "set title ", "'Bench Solution'"
      write(file_ref,*) "set xlabel 't'"
      write(file_ref,*) "set ylabel 'x'"
      ! write(file_ref,*) "set ylabel 'z'"
      write(file_ref,*) "set style line 1 lc rgb '#0060ad' lt 1 lw 1 pt 5 ps 1   # blue"
      write(file_ref,*) "plot '",input_file,"' using 1:2 with linespoints ls 1"
      close(file_ref)
      ! execute code
      call system("gnuplot " // output_file)
      !delete the plt file upon creation
      open(unit=1234, iostat=stat, file=output_file, status='old')
      if (stat == 0) close(1234, status='delete')
    else
      write(*,*) 'Error no data file exists'
    endif

  end subroutine


end module
