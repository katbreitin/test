module timer_mod


integer, parameter :: prec=kind(1d0), i64=selected_int_kind(15)
type Timer
   private

   integer :: start, rate=-1
   integer,allocatable :: class_start(:)
   real,allocatable :: class_rate(:)
   integer, allocatable :: class_count(:)
   character(len=200), allocatable:: class_list(:)
   integer :: n_class
   real, allocatable :: sec_per_class(:)
 contains
   procedure, public :: Tic, Tac, init, summary, reset
 end type Timer


contains

  subroutine  init(self, class_list)
    class (Timer), intent(inout) :: self
    character(len=*), intent(in) :: class_list(:)
    integer :: i

    self % n_class = size(class_list)

    allocate(self % sec_per_class(self % n_class))
    allocate(self % class_start(self % n_class))
    allocate(self % class_rate(self % n_class))
    allocate(self % class_count(self % n_class))
    allocate(self % class_list(self % n_class))
    self % class_list = class_list


  end subroutine

  subroutine Tic(self,class)
     class (Timer), intent(inout) :: self
     integer, intent(in) :: class
     integer :: start, rate

     call system_clock(count_rate=rate)
     call system_clock(start)
     self % start=start
     self % rate=rate

    self % class_start(class) = start
    self % class_rate(class) = rate
    self % class_count(class) = self % class_count(class) + 1
   end subroutine tic

   subroutine Tac(self, class)
      class (Timer), intent(inout) :: self
      integer, intent(in) :: class
      integer :: finish

      if(self%rate<0) then
        print*, 'Call to ''Tac'' subroutine must come after call to ''Tic'''
        stop
      endif

      call system_clock(finish)

      !print*, 'Elapsed time in seconds:', float(finish-self%start)/self%rate
      self % sec_per_class(class) =  self % sec_per_class(class) &
                    +  float(finish-self%class_start(class))/self%class_rate(class)



    end subroutine Tac


    subroutine summary(self, sort, minute)
       class (Timer), intent(inout) :: self
       logical, optional, intent(in) :: sort
       logical, optional, intent(in) :: minute
       real, allocatable:: sec_array(:)
       integer, allocatable:: sort_array(:)
       integer :: i
       integer :: idx
       character(len= 12) :: time_word
       character(len =50) :: FMT

        allocate(sec_array(self %n_class ))
        allocate(sort_array(self %n_class ))

       sec_array = self % sec_per_class

       if (present(sort)) then
           call quicksort(sec_array, sort_array)
       end if
       time_word = ' (sec): '
       if (present(minute)) then
          sec_array = sec_array/60.
          time_word = ' (minute): '
      end if
      FMT = "(I3,I3,I3, A, A, F10.5)"
      print*,'<--------   TIMING RESULTS ---------->'
      print*,'           RNK      NR      COUNT      DESCRIPTION                    TOTAL TIME '
       do i =  self % n_class , 1, -1
           idx = sort_array(i)

          write(*,FMT) self % n_class - i + 1,idx,self % class_count(idx),trim(self % class_list(idx)) &
             ,trim(time_word),sec_array(i)

       end do

       if (allocated(sec_array)) deallocate(sec_array)


    end subroutine

    subroutine reset (self)
       class (Timer), intent(inout) :: self

       self % sec_per_class(1:self%n_class) = 0
       self % class_rate(1:self%n_class) = 0
       self % class_start(1:self%n_class) = 0

    end  subroutine



    subroutine quicksort(array, sort_array)
    real, intent(inout)::array(:)
    real :: temp
    integer :: temp_sort
    integer :: i,j,last
    integer, intent(inout) :: sort_array(:)

    last=size(array)



    do i=1,last
      sort_array(i) = i
    end do


       do i=2,last
          temp=array(i)
          temp_sort = sort_array(i)
          do j=i-1,1,-1
             if (array(j).le.temp) exit
             array(j+1)=array(j)
             sort_array(j+1) = sort_array(j)
          enddo
          array(j+1)=temp
          sort_array(j+1) = temp_sort
       enddo



       return



  end subroutine quicksort



 end module timer_mod
