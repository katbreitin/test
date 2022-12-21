module cx_timer_mod
use timer_mod

type(timer) :: chronos_rttov


contains
subroutine timer_set_up (timer_obj)

   type(timer), intent(inout) :: timer_obj
   character(len=100), allocatable :: class_list(:)

   allocate (class_list(16))

   class_list(1)  =  "Level-1b Processing                 ."
   class_list(2)  =  "Ancil. Data Processing              ."
   class_list(3)  =  "RTM Processing                      ."
   class_list(4)  =  "Spatial Processing                  ."
   class_list(5)  =  "Default Aerosol                     ."
   class_list(6)  =  "Cloud Mask                          ."
   class_list(7)  =  "Cloud Type                          ."
   class_list(8)  =  "Cloud Height                        ."
   class_list(9)  =  "Lunar Cloud Opt/Micro               ."
   class_list(10) =  "Solar Cloud Opt/Micro               ."
   class_list(11) =  "Precip                              ."
   class_list(12) =  "Cloud Base and CCL                  ."
   class_list(13) =  "MURI Aerosol Retrieval              ."
   class_list(14) =  "Earth Radiation Budget              ."
   class_list(15) =  "Pixel-HDF Write                     ."
   class_list(16) =  "Total Time for Processing This Orbit."

   call timer_obj % init (class_list)

end subroutine


subroutine timer_set_up_all (timer_obj)

   type(timer), intent(inout) :: timer_obj
   character(len=100), allocatable :: class_list(:)

   allocate (class_list(1))
   class_list(1) =  "Total Time for All Orbits"
   call timer_obj % init (class_list)
end subroutine
end module
