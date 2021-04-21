#include "PersistentStream.h"

PersistentStream::PersistentStream(bool mode)
{
   m_mode = mode;
   m_stream = new strstream();
   m_frozen = false;
}

PersistentStream::~PersistentStream()
{
   if (m_frozen) m_stream->freeze(false);
   delete m_stream;
}

bool PersistentStream::IsSaving() const
{
   return m_mode;
}

void PersistentStream::write(const char *pch, int nCount)
{
   if (m_frozen) return;
   if (nCount > 0) m_stream->write(pch, nCount);
}

void PersistentStream::read(char *pch, int nCount)
{
   m_stream->read(pch, nCount);
}

void PersistentStream::clear()
{
   if (m_frozen)
   {
      m_stream->freeze(false);
      m_frozen = false;
   }
   delete m_stream;
   m_stream = new strstream();
}

int PersistentStream::size() const
{
   return m_stream->pcount();
}

const char *PersistentStream::GetBuf()
{
   m_frozen = true;
   return m_stream->str();
}

const bool PersistentStream::save = true;
const bool PersistentStream::load = false;
