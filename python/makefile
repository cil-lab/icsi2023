LIB_DIR = lib

default: icsi2023

icsi2023: setup.py icsi2023.pyx $(LIB_DIR)/icsi2023.a
	python setup.py build_ext --inplace && rm -f icsi2023.c && rm -Rf build

$(LIB_DIR)/icsi2023.a:
	make -C $(LIB_DIR) libicsi2023.a

clean:
	rm *.so