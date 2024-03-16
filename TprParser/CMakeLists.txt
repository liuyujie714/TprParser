cmake_minimum_required(VERSION 3.5)
project(TprParser)

# 设置编译选项
set(CMAKE_CXX_STANDARD 17)

# 包含 Python
find_package(PythonLibs 3 REQUIRED)
find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
if (WIN32)
	include_directories(${Python3_INCLUDE_DIRS} "${Python3_INCLUDE_DIRS}/../Lib/site-packages/numpy/core/include")
else()
	include_directories(${Python3_INCLUDE_DIRS} "${Python3_INCLUDE_DIRS}/../../lib/python3.8/site-packages/numpy/core/include")
endif()

# 添加可执行文件
add_executable(TprReader src/main.cpp src/Reader.cpp)

# 添加动态链接库
add_library(TprParser_ SHARED src/Reader.cpp src/Py_Module.cpp)

# windows or linux 
if(WIN32)
	set_target_properties(TprParser_ PROPERTIES SUFFIX ".pyd")
else()
	set_target_properties(TprParser_ PROPERTIES SUFFIX ".so")
	set_target_properties(TprParser_ PROPERTIES PREFIX "") # 去除lib前缀
endif()

# 链接 Python 库
target_link_libraries(TprParser_ ${PYTHON_LIBRARIES})

# 设置输出路径
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)

# 安装可执行文件和动态链接库
install(TARGETS TprReader DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
install(TARGETS TprParser_ DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)

