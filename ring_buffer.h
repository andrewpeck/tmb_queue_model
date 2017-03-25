#ifndef RING_BUFFER_H
#define RING_BUFFER_H 


class ringBuffer {

    public:

        ringBuffer (int ring_size)
        {
            _size = ring_size;
            buffer.reserve (ring_size);
        }

        void write (int data) {
            buffer[_wr_adr] = data;
        }

        int read  () {
            return buffer[_rd_adr];
        }

        void inc_wr_adr (int n_steps=1) {
            _wr_adr = (_wr_adr + n_steps) % _size;
        }

        void inc_rd_adr (int n_steps=1) {
            _rd_adr = (_rd_adr + n_steps) % _size;
        }

        int wr_adr () {
            return _wr_adr;
        }

        int rd_adr () {
            return _rd_adr;
        }

        void rd_adr (int adr) {
            _rd_adr = adr;
        }

    private:

        int _wr_adr=0;
        int _rd_adr=0;

        int _size=0;

        std::vector <int> buffer;

};

#endif /* RING_BUFFER_H */
