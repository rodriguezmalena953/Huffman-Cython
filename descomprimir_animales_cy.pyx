# cython: language_level=3

import heapq
import struct

ARCHIVO_ENTRADA = "animales.huf"
ARCHIVO_SALIDA = "animales_recuperado.txt"
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


cdef tuple leer_cabecera(bytes datos):
    cdef Py_ssize_t pos = 0
    cdef unsigned long long tamanio_original
    cdef unsigned int cantidad
    cdef unsigned int i
    cdef unsigned int simbolo
    cdef unsigned long long frecuencia
    cdef list frecuencias = [0] * MAX_SIMBOLOS
    cdef unsigned long long suma = 0

    if len(datos) < 14:
        raise ValueError("Archivo demasiado pequeno para ser HUF1")

    if datos[0:4] != b"HUF1":
        raise ValueError("El archivo no tiene formato HUF1")

    pos = 4
    tamanio_original = struct.unpack_from("<Q", datos, pos)[0]
    pos += 8

    cantidad = struct.unpack_from("<H", datos, pos)[0]
    pos += 2

    if cantidad > MAX_SIMBOLOS:
        raise ValueError("Cantidad de simbolos invalida")

    for i in range(cantidad):
        if pos + 9 > len(datos):
            raise ValueError("Cabecera incompleta")

        simbolo = datos[pos]
        pos += 1
        frecuencia = struct.unpack_from("<Q", datos, pos)[0]
        pos += 8

        if frecuencia == 0:
            raise ValueError("Frecuencia invalida")

        frecuencias[simbolo] = frecuencia

    for i in range(MAX_SIMBOLOS):
        suma += frecuencias[i]

    if suma != tamanio_original:
        raise ValueError("La suma de frecuencias no coincide con el tamano original")

    return tamanio_original, frecuencias, pos


cdef bytes descomprimir_datos(bytes datos_comprimidos, object raiz, unsigned long long tamanio_original):
    cdef bytearray salida = bytearray()
    cdef object actual
    cdef unsigned int byte
    cdef int bit
    cdef int valor_bit
    cdef Py_ssize_t i
    cdef const unsigned char[:] vista = datos_comprimidos

    if tamanio_original == 0:
        return b""

    if raiz is None:
        raise ValueError("Arbol Huffman invalido")

    if raiz.izq is None and raiz.der is None:
        return bytes([raiz.simbolo]) * tamanio_original

    actual = raiz

    for i in range(len(vista)):
        byte = vista[i]
        for bit in range(7, -1, -1):
            valor_bit = (byte >> bit) & 1
            actual = actual.der if valor_bit else actual.izq

            if actual is None:
                raise ValueError("Datos comprimidos corruptos")

            if actual.izq is None and actual.der is None:
                salida.append(actual.simbolo)

                if len(salida) == tamanio_original:
                    return bytes(salida)

                actual = raiz

    raise ValueError("Datos comprimidos incompletos")


def main():
    cdef bytes contenido
    cdef unsigned long long tamanio_original
    cdef list frecuencias
    cdef Py_ssize_t pos
    cdef object raiz
    cdef bytes recuperado

    try:
        with open(ARCHIVO_ENTRADA, "rb") as entrada:
            contenido = entrada.read()
    except OSError:
        print(f"Error: no se pudo abrir {ARCHIVO_ENTRADA}")
        return 1

    try:
        tamanio_original, frecuencias, pos = leer_cabecera(contenido)
        raiz = construir_arbol(frecuencias)
        recuperado = descomprimir_datos(contenido[pos:], raiz, tamanio_original)
    except ValueError as error:
        print(f"Error: {error}")
        return 1

    try:
        with open(ARCHIVO_SALIDA, "wb") as salida:
            salida.write(recuperado)
    except OSError:
        print(f"Error: no se pudo crear {ARCHIVO_SALIDA}")
        return 1

    print("Descompresion Huffman terminada.")
    print(f"Entrada : {ARCHIVO_ENTRADA}")
    print(f"Salida  : {ARCHIVO_SALIDA}")
    print(f"Bytes recuperados: {len(recuperado)}")

    return 0


if __name__ == "__main__":
    main()
