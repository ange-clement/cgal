// Copyright (c) 2020  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Maxime Gimeno

#ifdef CGAL_USE_SSH
#include "config.h"
#include <CGAL/Three/Three.h>

#include "CGAL/Use_ssh.h"
#include <CGAL/use.h>
#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <chrono>
#include <thread>
#include <sstream>

#include <QMessageBox>
#include <QStringList>


#include <libssh/sftp.h>
#include <fcntl.h>

bool test_result(int res)
{
  switch(res){

  case SSH_AUTH_ERROR:
  {
    std::cerr<<"Auth failed with error: "<<std::endl;
    return false;
  }
  case SSH_AUTH_DENIED:
  {
    std::cerr<<"The server doesn't accept that public key as an authentication token. Try another key or another method."<<std::endl;
    return false;
  }
  case SSH_AUTH_PARTIAL :
  {
    std::cerr<<"You've been partially authenticated, you still have to use another method."<<std::endl;
    return false;
  }
  case SSH_OK:
    return true;
  case SSH_EOF:
    std::cerr<<"key doesn't exist."<<std::endl;
    return false;
  default:
    return false;
  }
  return true;
}
namespace CGAL{
namespace ssh_internal{

bool establish_ssh_session(ssh_session &session,
                           const char* user,
                           const char* server,
                           const char* pub_key_path,
                           const char* priv_key_path,
                           const char* priv_key_password)
{
  int port = 22;

  //Can use SSH_LOG_PROTOCOL here for verbose output
  int verbosity = SSH_LOG_NOLOG;
  int res;
  //retry 4 times max each time the connection asks to be retried.
  for(int k = 0; k < 4; ++k)
  {
    if(session)
      ssh_free(session);
    session = ssh_new();
    ssh_options_set( session, SSH_OPTIONS_LOG_VERBOSITY, &verbosity );
    ssh_options_set( session, SSH_OPTIONS_PORT, &port );
    ssh_options_set( session, SSH_OPTIONS_USER, user );
    ssh_options_set( session, SSH_OPTIONS_HOST, server);

    ssh_connect(session);
#if LIBSSH_VERSION_MAJOR <1 && LIBSSH_VERSION_MINOR < 8
    if( ssh_is_server_known(session) != SSH_SERVER_KNOWN_OK )
#else
      if( ssh_session_is_known_server(session) != SSH_KNOWN_HOSTS_OK )
#endif
    {
      if(QMessageBox::warning(CGAL::Three::Three::mainWindow(), QString("Unknown Server"),
                              QString ("The server you are trying to join is not known.\n"
                                       "Do you wish to add it to the known servers list and continue?"),
                              QMessageBox::Yes | QMessageBox::No) != QMessageBox::Yes)
      {
        return false;
      }
#if LIBSSH_VERSION_MAJOR <1 && LIBSSH_VERSION_MINOR < 8
      if( ssh_write_knownhost(session) != SSH_OK )
#else
      if( ssh_session_update_known_hosts(session) != SSH_OK )
#endif
      {
        std::cerr << "writeKnownHost failed" << std::endl;
        return false;
      }
      else
      {
        ssh_connect(session);
      }
    }
    ssh_key pubkey = ssh_key_new();
    ssh_pki_import_pubkey_file(pub_key_path, &pubkey);
    res = ssh_userauth_try_publickey(session, nullptr, pubkey);
    ssh_key_free(pubkey);
    if(res == SSH_AUTH_AGAIN)
    {
      ssh_disconnect(session);
    }
    else
      break;
  }

  if(!test_result(res))
  {
    ssh_disconnect(session);

    return false;
  }

  ssh_key privkey = ssh_key_new();
  res = ssh_pki_import_privkey_file(priv_key_path, priv_key_password, nullptr, nullptr, &privkey);
  if (!test_result(res))
  {
    ssh_disconnect(session);
    ssh_key_free(privkey);
    return false;
  }
  res = ssh_userauth_publickey(session, nullptr, privkey);
  ssh_key_free(privkey);
  if(!test_result(res))
  {
    ssh_disconnect(session);
    return false;
  }
  return true;
}


bool establish_ssh_session_from_agent(ssh_session& session,
                                      const char *user,
                                      const char *server,
                                      const char *pub_key_path)
{
#ifndef _WIN32
  int port = 22;

  //Can use SSH_LOG_PROTOCOL here for verbose output
  int verbosity = SSH_LOG_NOLOG;
  int res;
  //retry 4 times max each time the connection asks to be retried.
  for(int k = 0; k < 4; ++k)
  {
    if(session)
      ssh_free(session);
    session = ssh_new();
    ssh_options_set( session, SSH_OPTIONS_LOG_VERBOSITY, &verbosity );
    ssh_options_set( session, SSH_OPTIONS_PORT, &port );
    ssh_options_set( session, SSH_OPTIONS_USER, user );
    ssh_options_set( session, SSH_OPTIONS_HOST, server);

    ssh_connect(session);
#if LIBSSH_VERSION_MAJOR <1 && LIBSSH_VERSION_MINOR < 8
    if( ssh_is_server_known(session) != SSH_SERVER_KNOWN_OK )
#else
    if( ssh_session_is_known_server(session) != SSH_KNOWN_HOSTS_OK )
#endif
    {
      if(QMessageBox::warning(CGAL::Three::Three::mainWindow(), QString("Unknown Server"),
                              QString ("The server you are trying to join is not known.\n"
                                       "Do you wish to add it to the known servers list and continue?"),
                              QMessageBox::Yes | QMessageBox::No) != QMessageBox::Yes)
      {
        return false;
      }
#if LIBSSH_VERSION_MAJOR <1 && LIBSSH_VERSION_MINOR < 8
      if( ssh_write_knownhost(session) != SSH_OK )
#else
      if( ssh_session_update_known_hosts(session) != SSH_OK )
#endif
      {
        std::cerr << "writeKnownHost failed" << std::endl;
        return false;
      }
      else
      {
        ssh_connect(session);
      }
    }
    ssh_key pubkey = ssh_key_new();
    ssh_pki_import_pubkey_file(pub_key_path, &pubkey);
    res = ssh_userauth_try_publickey(session, nullptr, pubkey);
    ssh_key_free(pubkey);
    if(res == SSH_AUTH_AGAIN)
      ssh_disconnect(session);
    else
      break;
  }


  if(!test_result(res))
  {
    ssh_disconnect(session);
    return false;
  }

  res = ssh_userauth_agent(session, user);
  if(!test_result(res))
  {
    ssh_disconnect(session);
    return false;
  }
  return true;
#else
  CGAL_USE(session);
  CGAL_USE(user);
  CGAL_USE(server);
  CGAL_USE(pub_key_path);

  return false;
#endif
}

void close_connection(ssh_session &session)
{
  ssh_disconnect(session);
}

bool push_file(ssh_session &session,
               const char* dest_path,
               const char* filepath)
{
  std::size_t processed = 0;
  sftp_file sftpfile;
  //copy a file
  sftp_session sftp = sftp_new(session);
  if (sftp == nullptr)
  {
    std::cerr<<"Error allocating sftp session:\n"
            << ssh_get_error(session)<<std::endl;
    return false;
  }
  int res = sftp_init(sftp);
  if(res < 0)
  {
    std::cerr<< "Error initializing sftp session:\n"
             << ssh_get_error(session)<<std::endl;
    sftp_free(sftp);
    ssh_disconnect(session);
    return false;
  }
  //read a file into a buffer
  std::ifstream file(filepath, std::ios::binary | std::ios::ate);
  if(!file)
  {
    std::cerr<<"File not found."<<std::endl;
    sftp_free(sftp);
    ssh_disconnect(session);
    return false;
  }
  std::streamsize size = file.tellg();
  file.seekg(0, std::ios::beg);

  std::vector<char> buffer(size);
  if (!file.read(buffer.data(), size))
  {
    std::cerr<<"error while reading file."<<std::endl;
    sftp_free(sftp);
    ssh_disconnect(session);
    return false;
  }
  //push a file to /tmp
  sftpfile = sftp_open(sftp, dest_path, O_WRONLY | O_CREAT, 0644);
  if (sftpfile == NULL)
  {
    std::cerr<< "Can't open remote file:\n"
             << ssh_get_error(session)<<std::endl;
    sftp_free(sftp);
    ssh_disconnect(session);
    return false;
  }
  while ( size > 0 )
  {
    int s = size;
    if (s > 16384)
      s = 16384;
    res = sftp_write(sftpfile, buffer.data() + processed, s);
    if ( res < 0)
    {
      std::cerr<< "Can't write data to file:\n"
               << ssh_get_error(session)<<std::endl;
      sftp_free(sftp);
      ssh_disconnect(session);
      return false;
    }
    size -= res;
    processed += res;
  }
  sftp_close(sftpfile);
  sftp_free(sftp);
  return true;
}

bool pull_file(ssh_session &session,
               const char* from_path,
               const char* to_path)
{
  std::size_t size;
  std::size_t processed = 0;
  std::vector<char> buffer;
  sftp_file sftpfile;
  sftp_session sftp = sftp_new(session);
  if (sftp == nullptr)
  {
    std::cerr<<"Error allocating sftp session:\n"
            << ssh_get_error(session)<<std::endl;
    return false;
  }
  int res = sftp_init(sftp);
  if(res < 0)
  {
    std::cerr<< "Error initializing sftp session:\n"
             << ssh_get_error(session)<<std::endl;
    sftp_free(sftp);
    ssh_disconnect(session);
    return false;
  }
  sftpfile = sftp_open(sftp, from_path, O_RDONLY, 0);
  if (sftpfile == NULL)
  {
    std::cerr<< "Can't open remote file:\n"
             << ssh_get_error(session)<<std::endl;
    sftp_free(sftp);
    ssh_disconnect(session);
    return false;
  }
  sftp_attributes sftpattr;
  sftpattr = sftp_stat(sftp,from_path);
  size=sftpattr->size;
  buffer.resize(size);
  while ( size > 0 )
  {
    int s = size;
    if (s > 16384)
      s = 16384;
    res = sftp_read(sftpfile, buffer.data() + processed, s);
    if ( res < 0)
    {
      std::cerr<< "Can't read data to file:\n"
               << ssh_get_error(session)<<std::endl;
      sftp_free(sftp);
      ssh_disconnect(session);
      return false;
    }
    size -= res;
    processed += res;
  }
  size=sftpattr->size;
  std::ofstream file(to_path, std::ios::binary |std::ios::trunc);
  if(!file.write(buffer.data(), size))
  {
    std::cerr<<"Error while writing file."<<std::endl;
  }
  file.close();
  sftp_close(sftpfile);
  sftp_free(sftp);
  return true;
}

bool explore_the_galaxy(ssh_session &session,
                        QStringList& files)
{
  ssh_channel channel;
  channel = ssh_channel_new(session);
  if (channel == nullptr) return false;
  int rc = ssh_channel_open_session(channel);
  if (rc != SSH_OK)
  {
    ssh_channel_free(channel);
    return rc;
  }
  rc = ssh_channel_request_exec(channel, "ls /tmp");
  if (rc != SSH_OK)
  {
    ssh_channel_close(channel);
    ssh_channel_free(channel);
    return rc;
  }

  char buffer[256];
  int nbytes;
  nbytes = ssh_channel_read(channel, buffer, sizeof(buffer), 0);
  while (nbytes > 0)
  {

    std::string sbuf(buffer, nbytes);
    if(sbuf.find("CGAL_Lab_") != std::string::npos)
    {
      std::istringstream iss(sbuf);
      std::string file;
      while(iss >> file)
      {
        if(file.find("CGAL_Lab_") != std::string::npos)
        {
          QString name(file.c_str());
          files.push_back(name.remove("CGAL_Lab_"));
        }
      }
    }

    nbytes = ssh_channel_read(channel, buffer, sizeof(buffer), 0);
  }
  if (nbytes < 0)
  {
    ssh_channel_close(channel);
    ssh_channel_free(channel);
    return false;
  }
  ssh_channel_send_eof(channel);
  ssh_channel_close(channel);
  ssh_channel_free(channel);
  return true;
}

}// end of ssh_internal
}// end of CGAL
#endif
