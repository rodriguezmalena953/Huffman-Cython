# cython: language_level=3

import heapq
import os
import struct

ARCHIVO_ENTRADA = "animales.txt"
ARCHIVO_SALIDA = "animales.huf"
MAX_SIMBOLOS = 256


cdef class NodoHuffman:
    cdef public int simbolo
    cdef public unsigned long long frecuencia
    cdef public unsigned int orden
    cdef public object izq
    cdef public object der

    def __cinit__(self, int simbolo, unsigned long long frecuencia, unsigned int orden):
        self.simbolo = simbolo
        self.frecuencia = frecuencia
        self.orden = orden
        self.izq = None
        self.der = None


cdef object construir_arbol(list frecuencias):
    cdef list heap = []
    cdef unsigned int orden = 0
    cdef int i
    cdef object a, b, padre

    for i in range(MAX_SIMBOLOS):
        if frecuencias[i] > 0:
            nodo = NodoHuffman(i, frecuencias[i], orden)
            heapq.heappush(heap, (frecuencias[i], orden, nodo))
            orden += 1

    if not heap:
        return None

    while len(heap) > 1:
        _, _, a = heapq.heappop(heap)
        _, _, b = heapq.heappop(heap)

        padre = NodoHuffman(-1, a.frecuencia + b.frecuencia, orden)
        orden += 1
        padre.izq = a
        padre.der = b

        heapq.heappush(heap, (padre.frecuencia, padre.orden, padre))

    return heapq.heappop(heap)[2]


cdef void generar_codigos(object nodo, str camino, dict codigos):
    if nodo is None:
        return

    if nodo.izq is None and nodo.der is None:
        if camino == "":
            camino = "0"
        codigos[nodo.simbolo] = camino
        return

    generar_codigos(nodo.izq, camino + "0", codigos)
    generar_codigos(nodo.der, camino + "1", codigos)


cdef bytes crear_cabecera(unsigned long long tamanio_original, list frecuencias):
    cdef bytearray cabecera = bytearray()
    cdef int i
    cdef unsigned int cantidad = 0

    cabecera.extend(b"HUF1")
    cabecera.extend(struct.pack("<Q", tamanio_original))

    for i in range(MAX_SIMBOLOS):
        if frecuencias[i] > 0:
            cantidad += 1

    cabecera.extend(struct.pack("<H", cantidad))

    for i in range(MAX_SIMBOLOS):
        if frecuencias[i] > 0:
            cabecera.append(i)
            cabecera.extend(struct.pack("<Q", frecuencias[i]))

    return bytes(cabecera)


cdef bytes comprimir_datos(bytes datos, dict codigos):
    cdef bytearray salida = bytearray()
    cdef unsigned int byte_salida = 0
    cdef int bits_usados = 0
    cdef unsigned int simbolo
    cdef Py_ssize_t i
    cdef str codigo
    cdef str bit
    cdef const unsigned char[:] vista = datos

    for i in range(len(vista)):
        simbolo = vista[i]
        codigo = codigos[simbolo]

        for bit in codigo:
            byte_salida = (byte_salida << 1) | (1 if bit == "1" else 0)
            bits_usados += 1

            if bits_usados == 8:
                salida.append(byte_salida)
                byte_salida = 0
                bits_usados = 0

    if bits_usados > 0:
        byte_salida <<= (8 - bits_usados)
        salida.append(byte_salida)

    return bytes(salida)


def main():
    cdef list frecuencias = [0] * MAX_SIMBOLOS
    cdef bytes datos
    cdef unsigned int simbolo
    cdef Py_ssize_t i
    cdef const unsigned char[:] vista
    cdef unsigned long long tamanio_original
    cdef object raiz
    cdef dict codigos = {}
    cdef bytes cabecera
    cdef bytes comprimidos
    cdef long tamanio_comprimido

    try:
        with open(ARCHIVO_ENTRADA, "rb") as entrada:
            datos = entrada.read()
    except OSError:
        print(f"Error: no se pudo abrir {ARCHIVO_ENTRADA}")
        return 1

    tamanio_original = len(datos)

    vista = datos
    for i in range(len(vista)):
        simbolo = vista[i]
        frecuencias[simbolo] += 1

    raiz = construir_arbol(frecuencias)
    if raiz is not None:
        generar_codigos(raiz, "", codigos)

    cabecera = crear_cabecera(tamanio_original, frecuencias)
    comprimidos = comprimir_datos(datos, codigos) if tamanio_original > 0 else b""

    try:
        with open(ARCHIVO_SALIDA, "wb") as salida:
            salida.write(cabecera)
            salida.write(comprimidos)
    except OSError:
        print(f"Error: no se pudo crear {ARCHIVO_SALIDA}")
        return 1

    tamanio_comprimido = os.path.getsize(ARCHIVO_SALIDA)

    print("Compresion Huffman terminada.")
    print(f"Entrada : {ARCHIVO_ENTRADA}")
    print(f"Salida  : {ARCHIVO_SALIDA}")
    print(f"Tamano original   : {tamanio_original} bytes")
    print(f"Tamano comprimido : {tamanio_comprimido} bytes")

    if tamanio_original > 0:
        print(f"Tamano final      : {100.0 * tamanio_comprimido / tamanio_original:.2f} % del original")

    return 0


if __name__ == "__main__":
    main()
