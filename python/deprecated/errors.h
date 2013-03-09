#define PY_CATCHALL \
    catch (tmv::Error& e) { \
      throw e.what(); \
    } catch (FileNotFoundException& e) { \
      throw e.what(); \
    } catch (ParameterException& e) { \
      throw e.what(); \
    } catch (ReadException& e) { \
      throw e.what(); \
    } catch (WriteException& e) { \
      throw e.what(); \
    } catch (ProcessingException& e) { \
      throw e.what(); \
    } catch (std::exception& e) { \
      throw e.what(); \
    } catch (...) { \
        throw "Caught unknown error"; \
    } 

